"""ESO/MPG 2.2 m Wide Field Imager: calibrate, plate solve and extract one night.

Usage
-----
    python -m pipeline.wfi calibrate SEQ_DIR
    python -m pipeline.wfi solve SEQ_DIR
    python -m pipeline.wfi extract SEQ_DIR --chip N --planet "POTS-1 b" --out OUT.json [--plot OUT.png]

``SEQ_DIR`` is a directory written by ``python -m access.fetch`` (``science/``,
``calib/``, ``manifest.json``). Products go in ``SEQ_DIR/products/``.

WFI raw frames are 8 image extensions, one per CCD (2046 x 4098 imaging
pixels). Each has 48 columns of prescan and 48 of overscan in x, and 30
overscan rows in y. The y rows are at the bottom of some chips and the top of
others, so their side is found from the data. The header gives gain for two
chips only and read noise for none, so both are measured here from the
night's biases and screen flats (photon transfer). Raw frames have no
celestial WCS. Each chip of the first science frame is solved against Gaia.
"""

import argparse
import json
import logging
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

from exotic.extract.apertures import FwhmScaled
from exotic.extract.sequence import Observatory, extract_sequence

from . import store
from .astrometry import gaia_stars, solve_wcs_robust
from .lightcurve import LIMB_DARKENING, choose_sources, fit_with_airmass, pick_comparisons

log = logging.getLogger("pipeline.wfi")

PIXEL_SCALE = 0.238          # arcsec/px, nominal
OVERSCAN_ROWS = 30
SATURATION = 50000.0         # ADU after calibration; raw frames clip at 65535
N_CHIPS = 8
# Solve attempts for one 8' x 16' chip: (Gaia circle radius deg, detections, Gaia stars).
CHIP_SOLVE_ATTEMPTS = ((0.12, 80, 100), (0.12, 50, 60), (0.15, 80, 140), (0.10, 60, 60))


# --- geometry and per-frame calibration ------------------------------------------------

def _det_value(header, name):
    """An ``ESO DET OUTn <name>`` value, whatever this extension's output number n is."""
    for key in header:
        if re.fullmatch(rf"ESO DET OUT\d {name}", key):
            return header[key]
    raise KeyError(name)


def chip_number(header):
    return int(re.search(r"CHIP(\d)", header["EXTNAME"]).group(1))


def overscan_corrected(data, header):
    """Overscan-subtracted imaging section of one raw chip, as float32.

    Subtracts the per-row median of the x overscan, then drops the prescan,
    overscan and the 30 y-overscan rows. The y rows are whichever end is
    closer to the x-overscan level.
    """
    data = np.asarray(data, dtype=np.float32)
    prscx, nx, ovscx = (int(_det_value(header, k)) for k in ("PRSCX", "NX", "OVSCX"))
    overscan = data[:, prscx + nx + 4:prscx + nx + ovscx - 4]
    level = np.median(overscan, axis=1)
    image = data[:, prscx:prscx + nx] - level[:, None]
    bottom = np.median(image[:OVERSCAN_ROWS])
    top = np.median(image[-OVERSCAN_ROWS:])
    return image[OVERSCAN_ROWS:] if abs(bottom) < abs(top) else image[:-OVERSCAN_ROWS]


def read_chip(path, ext):
    with fits.open(path) as hdul:  # BZERO-scaled ints, so no memmap
        return overscan_corrected(hdul[ext].data, hdul[ext].header), hdul[ext].header.copy(), hdul[0].header.copy()


# --- calibration products ----------------------------------------------------------------

def _calib_files(seq_dir, dp_type):
    manifest = json.loads((seq_dir / "manifest.json").read_text())
    entries = [c for c in manifest["calibrations"] if c["dp_type"] == dp_type]
    return [(seq_dir / "calib" / f"{c['dp_id']}.fits", c) for c in entries
            if (seq_dir / "calib" / f"{c['dp_id']}.fits").exists()]


def _centre(image, half=256):
    ny, nx = image.shape
    return image[ny // 2 - half:ny // 2 + half, nx // 2 - half:nx // 2 + half]


def calibrate(seq_dir):
    """Master bias and flat per chip, and measured gain / read noise. Writes ``products/``."""
    seq_dir = Path(seq_dir)
    out = seq_dir / "products"
    out.mkdir(exist_ok=True)
    biases = [p for p, _ in _calib_files(seq_dir, "BIAS")]
    sky_flats = [(p, c) for p, c in _calib_files(seq_dir, "FLAT,SKY")]
    screen = [(p, c) for p, c in _calib_files(seq_dir, "FLAT,SCREEN")]
    log.info("%d biases, %d sky flats, %d screen flats", len(biases), len(sky_flats), len(screen))

    summary = {}
    for ext in range(1, N_CHIPS + 1):
        bias_frames = np.array([read_chip(p, ext)[0] for p in biases])
        header = fits.getheader(biases[0], ext)
        chip = chip_number(header)
        master_bias = np.median(bias_frames, axis=0)

        # Read noise from bias pairs; gain from same-exposure screen-flat pairs (photon transfer).
        ron_adu = np.median([sigma_clipped_stats(_centre(a - b))[2] / np.sqrt(2)
                             for a, b in zip(bias_frames[0::2], bias_frames[1::2])])
        by_exptime = {}
        for p, c in screen:
            by_exptime.setdefault(round(c["exposure"], 1), []).append(read_chip(p, ext)[0] - master_bias)
        gains = []
        for frames in by_exptime.values():
            for a, b in zip(frames[0::2], frames[1::2]):
                signal = np.median(_centre(a)) + np.median(_centre(b))
                var = sigma_clipped_stats(_centre(a - b))[2] ** 2 - 2 * ron_adu ** 2
                gains.append(signal / var)
        gain = float(np.median(gains))

        # Flat: normalised median of the twilight sky flats in their linear range, else the screen flats.
        usable = []
        for p, _ in sky_flats:
            f = read_chip(p, ext)[0] - master_bias
            level = np.median(_centre(f))
            if 3000 < level < 35000:
                usable.append(f / level)
        source = "sky"
        if len(usable) < 3:
            usable = [f / np.median(_centre(f)) for frames in by_exptime.values() for f in frames]
            source = "screen"
        master_flat = np.median(usable, axis=0)
        master_flat /= np.median(master_flat)
        master_flat[master_flat < 0.2] = np.nan  # dead / vignetted pixels

        fits.writeto(out / f"chip{chip}_bias.fits", master_bias.astype(np.float32), overwrite=True)
        fits.writeto(out / f"chip{chip}_flat.fits", master_flat.astype(np.float32), overwrite=True)
        summary[chip] = {"ext": ext, "gain_e_per_adu": gain, "read_noise_e": float(ron_adu * gain),
                         "read_noise_adu": float(ron_adu), "n_gain_pairs": len(gains),
                         "gain_pair_scatter": float(np.std(gains)), "flat_source": source,
                         "n_flats": len(usable), "bias_level_adu": float(np.median(master_bias)),
                         "header_gain": header.get("HIERARCH ESO DET OUT1 CONAD")}
        log.info("chip %d (ext %d): gain %.3f e/ADU (%d pairs, scatter %.3f), RON %.2f e, flat from %d %s",
                 chip, ext, gain, len(gains), np.std(gains), ron_adu * gain, len(usable), source)
    (out / "detector.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def calibrated_chip(path, ext, bias, flat):
    image, chip_header, primary = read_chip(path, ext)
    return (image - bias) / flat, chip_header, primary


def _load_products(seq_dir, chip):
    out = Path(seq_dir) / "products"
    detector = json.loads((out / "detector.json").read_text())[str(chip)]
    return (detector, fits.getdata(out / f"chip{chip}_bias.fits"), fits.getdata(out / f"chip{chip}_flat.fits"))


def science_files(seq_dir):
    manifest = json.loads((Path(seq_dir) / "manifest.json").read_text())
    files = [Path(seq_dir) / "science" / f"{s['dp_id']}.fits" for s in manifest["science"]]
    return [p for p in files if p.exists()]


def site(primary):
    return Observatory(latitude=float(primary["HIERARCH ESO TEL GEOLAT"]),
                       longitude=float(primary["HIERARCH ESO TEL GEOLON"]),
                       elevation=float(primary["HIERARCH ESO TEL GEOELEV"]))


# --- astrometry --------------------------------------------------------------------------

def chip_centre_guesses(primary, chip_header):
    """Guessed sky centres of a chip, from the pointing and the chip's place in the 4 x 2 mosaic.

    The mosaic's orientation on the sky isn't in the header, so all four
    sign combinations are tried, then the pointing itself.
    """
    chip = chip_number(chip_header)
    cx = chip_header.get(f"HIERARCH ESO DET CHIP{chip} X")
    cy = chip_header.get(f"HIERARCH ESO DET CHIP{chip} Y")
    ra0, dec0 = float(primary["RA"]), float(primary["DEC"])
    guesses = []
    if cx is not None and cy is not None:
        dx = (cx - 2.5) * 2046 * PIXEL_SCALE / 3600
        dy = (cy - 1.5) * 4098 * PIXEL_SCALE / 3600
        for sx in (1, -1):
            for sy in (1, -1):
                for ddx, ddy in ((sx * dx, sy * dy), (sy * dy, sx * dx)):  # with and without a 90 deg turn
                    guesses.append((ra0 + ddx / np.cos(np.radians(dec0)), dec0 + ddy))
    guesses.append((ra0, dec0))
    return list(dict.fromkeys((round(r, 4), round(d, 4)) for r, d in guesses))


def solve(seq_dir, g_max=17.0):
    """WCS for every chip of the first science frame. Writes ``products/chipN_wcs.fits`` and ``wcs.json``."""
    seq_dir = Path(seq_dir)
    out = seq_dir / "products"
    reference = science_files(seq_dir)[0]
    primary = fits.getheader(reference, 0)
    catalogue = gaia_stars(float(primary["RA"]), float(primary["DEC"]), 0.45, g_max)
    log.info("Gaia catalogue for the mosaic: %d stars G < %.1f", len(catalogue), g_max)

    summary = {}
    for ext in range(1, N_CHIPS + 1):
        chip_header = fits.getheader(reference, ext)
        chip = chip_number(chip_header)
        _, bias, flat = _load_products(seq_dir, chip)
        image = np.nan_to_num(calibrated_chip(reference, ext, bias, flat)[0])
        try:
            wcs, _, stats = solve_wcs_robust(image, chip_centre_guesses(primary, chip_header), PIXEL_SCALE,
                                             g_max=g_max, fwhm=4.0, match_px=4.0, min_matched=30,
                                             max_rms_arcsec=0.5, attempts=CHIP_SOLVE_ATTEMPTS, catalogue=catalogue)
        except RuntimeError as exc:
            log.warning("chip %d: %s", chip, str(exc).splitlines()[0])
            summary[chip] = {"ext": ext, "solved": False}
            continue
        fits.PrimaryHDU(header=wcs.to_header(relax=True)).writeto(out / f"chip{chip}_wcs.fits", overwrite=True)
        centre = wcs.pixel_to_world(image.shape[1] / 2, image.shape[0] / 2)
        summary[chip] = {"ext": ext, "solved": True, **stats, "centre": [centre.ra.deg, centre.dec.deg]}
        log.info("chip %d: %d stars, rms %.3f arcsec, centre %.4f %.4f", chip, stats["n_matched"],
                 stats["rms_arcsec"], centre.ra.deg, centre.dec.deg)
    (out / "wcs.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    return summary


def chip_containing(seq_dir, coord, margin_px=30):
    for path in sorted((Path(seq_dir) / "products").glob("chip*_wcs.fits")):
        wcs = WCS(fits.getheader(path))
        x, y = wcs.world_to_pixel(coord)
        if margin_px < x < 2046 - margin_px and margin_px < y < 4098 - OVERSCAN_ROWS - margin_px:
            return int(re.search(r"chip(\d)", path.name).group(1)), float(x), float(y)
    return None


# --- extraction and fit ------------------------------------------------------------------

def planet_parameters(name, bjd_mid):
    """Orbit and ephemeris from the NASA Exoplanet Archive, projected to this night."""
    import pyvo

    row = pyvo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP").search(
        "SELECT pl_name, ra, dec, pl_orbper, pl_orbpererr1, pl_tranmid, pl_tranmiderr1, pl_ratror, pl_ratdor, "
        "pl_orbsmax, st_rad, pl_orbincl, pl_orbeccen, pl_orblper, sy_gaiamag "
        f"FROM pscomppars WHERE pl_name='{name}'").to_table()[0]

    def val(key, default=0.0):
        v = row[key]
        return default if np.ma.is_masked(v) or not np.isfinite(v) else float(v)

    per, t0 = val("pl_orbper"), val("pl_tranmid")
    ars = val("pl_ratdor", np.nan)
    if not np.isfinite(ars):  # not always tabulated: a / R_star from the semi-major axis and stellar radius
        ars = val("pl_orbsmax") * 215.032 / val("st_rad", 1.0)
    epoch = int(np.round((bjd_mid - t0) / per))
    sigma = np.hypot(val("pl_tranmiderr1"), abs(epoch) * val("pl_orbpererr1"))
    return {"name": name, "ra": val("ra"), "dec": val("dec"), "per": per, "ars": ars,
            "inc": val("pl_orbincl", 90.0), "ecc": val("pl_orbeccen"), "omega": val("pl_orblper"),
            "rprs": val("pl_ratror", 0.1), "t0_published": t0, "epoch": epoch, "tmid_predicted": t0 + epoch * per,
            "tmid_predicted_sigma_min": float(sigma * 1440), "source": "NASA Exoplanet Archive pscomppars"}


def _for_exotic(image):
    """A calibrated frame in the form EXOTIC's background estimate needs: whole ADU.

    ``exotic.extract.photometry.skybg_phot`` takes the statistical mode of the
    pixel values around each star. On flat-fielded float data almost every
    value is unique, so the mode is an arbitrary pixel and the background jumps
    by tens of ADU from frame to frame (20-30% flux errors at G~17 on WFI).
    Rounding to integers restores a meaningful mode, at a cost of 0.29 ADU rms
    of quantisation noise against ~7 ADU of sky noise.
    """
    return np.round(np.nan_to_num(image, nan=0.0))


def _extract_chunk(job):
    """``extract_sequence`` for a subset of sources, in a worker process."""
    detector, bias, flat = _load_products(job["seq_dir"], job["chip"])

    def frames():
        for i, path in enumerate(job["files"]):
            image, _, header = calibrated_chip(path, detector["ext"], bias, flat)
            if i == 0:
                header.update(job["wcs_header"])
            yield _for_exotic(image), header

    return extract_sequence(frames(), sky_positions=job["sky_positions"], ra=job["ra"], dec=job["dec"],
                            observatory=job["observatory"], apertures=job["apertures"], gain=job["gain"],
                            read_noise=job["read_noise"])


def extract(seq_dir, chip, planet_name, out_path, plot_path=None, g_max=18.5, max_sources=800,
            aperture_coeff=None, workers=16):
    seq_dir = Path(seq_dir)
    detector, bias, flat = _load_products(seq_dir, chip)
    ext = detector["ext"]
    wcs = WCS(fits.getheader(seq_dir / "products" / f"chip{chip}_wcs.fits"))
    files = science_files(seq_dir)
    reference, _, primary = calibrated_chip(files[0], ext, bias, flat)
    obs = site(primary)

    t_mid_guess = float(primary["MJD-OBS"]) + 2400000.5
    planet = planet_parameters(planet_name, t_mid_guess)
    target = SkyCoord(planet["ra"] * u.deg, planet["dec"] * u.deg)

    ny, nx = reference.shape
    centre = wcs.pixel_to_world(nx / 2, ny / 2)
    gaia = gaia_stars(centre.ra.deg, centre.dec.deg, np.hypot(nx, ny) / 2 * PIXEL_SCALE / 3600 * 1.02, g_max + 2.5)
    picked = choose_sources(np.nan_to_num(reference), wcs, gaia, target, g_max, edge_px=30,
                            blend_px=4.0 / PIXEL_SCALE, saturation=SATURATION, max_sources=max_sources)
    gaia = gaia[picked]
    log.info("chip %d: photometering %d sources (target G=%.2f)", chip, len(gaia), gaia["phot_g_mean_mag"][0])

    # Radii are fixed once here, on the whole source list, so every worker uses the same ones.
    apertures = FwhmScaled() if aperture_coeff is None else FwhmScaled(aperture_coeff=aperture_coeff)
    ref_x, ref_y = wcs.world_to_pixel(SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg))
    radii = apertures.resolve(_for_exotic(reference), list(zip(ref_x, ref_y)))
    log.info("apertures: r=%.2f, annulus %.2f-%.2f px", radii.aperture_r, radii.annulus_r_in, radii.annulus_r_out)

    # Centroid fitting dominates (~70 ms per source per frame), so split the
    # sources across processes. Each one re-registers the frames itself.
    sky_positions = list(zip(np.asarray(gaia["ra"], float), np.asarray(gaia["dec"], float)))
    chunks = [sky_positions[i::workers] for i in range(workers)]
    order = np.concatenate([np.arange(len(sky_positions))[i::workers] for i in range(workers)])
    job = dict(seq_dir=str(seq_dir), chip=chip, files=[str(f) for f in files], wcs_header=wcs.to_header(relax=True),
               ra=target.ra.deg, dec=target.dec.deg, observatory=obs, apertures=radii,
               gain=detector["gain_e_per_adu"], read_noise=detector["read_noise_e"])
    with ProcessPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(_extract_chunk, [dict(job, sky_positions=c) for c in chunks if c]))
    grid = np.concatenate([p.reshape(len(files), -1) for p in parts], axis=1)[:, np.argsort(order)]
    flux = grid["flux"]
    bjd, airmass = grid["bjd_tdb"][:, 0], grid["airmass"][:, 0]
    aligned = grid["aligned"].all(axis=1)

    good = aligned & np.isfinite(flux[:, 0]) & (flux[:, 0] > 0)
    comps, comp_scatter = pick_comparisons(flux[good][:, 1:], flux[good][:, 0], n_max=15, ratio_range=(0.5, 5.0))
    comps = [c + 1 for c in comps]
    relative = flux[:, 0] / flux[:, comps].sum(axis=1)
    relative = relative / np.nanmedian(relative[good])
    values, errors, extra = fit_with_airmass(bjd[good], relative[good], airmass[good], planet)
    fit = {"tmid_bjd": values["tmid"], "tmid_err_min": errors["tmid"] * 1440, "rprs": values["rprs"],
           "rprs_err": errors["rprs"], "scale": values["scale"],
           "o_minus_c_min": (values["tmid"] - planet["tmid_predicted"]) * 1440,
           "airmass_poly": extra["airmass_poly"], "residual_rms": extra["rms"]}
    log.info("fit: T_mid %.5f +- %.2f min, O-C %+.2f min, Rp/Rs %.3f +- %.3f, rms %.4f", fit["tmid_bjd"],
             fit["tmid_err_min"], fit["o_minus_c_min"], fit["rprs"], fit["rprs_err"], fit["residual_rms"])

    roles = ["target"] + ["comparison" if j in comps else "field" for j in range(1, len(gaia))]
    sources = [{"gaia_dr3": int(gaia["source_id"][j]), "g_mag": float(gaia["phot_g_mean_mag"][j]),
                "ra": float(gaia["ra"][j]), "dec": float(gaia["dec"][j]), "role": roles[j],
                "flux": grid["flux"][:, j], "uncertainty": grid["uncertainty"][:, j],
                "x": grid["x_centroid"][:, j], "y": grid["y_centroid"][:, j]} for j in range(len(gaia))]
    manifest = json.loads((seq_dir / "manifest.json").read_text())
    sequence = {"archive": "ESO", "instrument": "WFI", "chip": chip, "object": primary.get("OBJECT"),
                "filter": primary.get("HIERARCH ESO INS FILT1 NAME"), "exposure_s": float(primary["EXPTIME"]),
                "prog_id": manifest["target"].get("prog_id"), "night": manifest["target"].get("night"),
                "dp_ids": [p.stem for p in files],
                "site": {"latitude": obs.latitude, "longitude": obs.longitude, "elevation": obs.elevation},
                "apertures": {k: float(grid[k][0, 0]) for k in ("aperture_r", "annulus_r_in", "annulus_r_out")},
                "gain": detector["gain_e_per_adu"], "read_noise": detector["read_noise_e"],
                "uncertainties": "measured gain/read noise (photon transfer)", "detector": detector,
                "astrometry": json.loads((seq_dir / "products" / "wcs.json").read_text())[str(chip)],
                "g_max": g_max, "limb_darkening": LIMB_DARKENING}
    frames_out = {"bjd_tdb": bjd, "jd_utc": grid["jd_utc"][:, 0], "airmass": airmass, "aligned": aligned,
                  "used": good}
    analysis = {"planet": planet, "comparisons": {"gaia_dr3": [int(gaia["source_id"][c]) for c in comps],
                                                  "point_to_point_scatter": list(comp_scatter.values())},
                "relative_flux": relative, "transit_fit": fit}
    store.write(out_path, sequence, frames_out, sources, analysis)
    log.info("wrote %s", out_path)

    if plot_path:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = (bjd[good] - planet["tmid_predicted"]) * 24
        order = np.argsort(t)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(t, extra["detrended"], ".", color="0.4", ms=4, label="target / comparison ensemble, detrended")
        nbin = 5
        nb = len(t) // nbin * nbin
        ax.plot(t[order][:nb].reshape(-1, nbin).mean(1), extra["detrended"][order][:nb].reshape(-1, nbin).mean(1),
                "o", color="C0", ms=5, label=f"binned x{nbin}")
        ax.plot(t[order], extra["model"][order], "-", color="C3", label="fit_transit")
        ax.axvline(0, color="0.6", ls="--", lw=1, label="predicted T_mid")
        ax.set_xlabel("hours from predicted mid-transit")
        ax.set_ylabel("relative flux")
        ax.set_title(f"{planet['name']} (WFI chip {chip}): T_mid {fit['tmid_bjd']:.5f} ± {fit['tmid_err_min']:.1f} min, "
                     f"O-C {fit['o_minus_c_min']:+.1f} min, Rp/Rs {fit['rprs']:.3f}", fontsize=10)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=120)
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("calibrate", "solve", "extract"):
        p = sub.add_parser(name)
        p.add_argument("seq_dir")
        if name == "extract":
            p.add_argument("--chip", type=int, help="chip number; default: the chip containing the planet's host")
            p.add_argument("--planet", required=True, help="NASA Exoplanet Archive pl_name, e.g. 'POTS-1 b'")
            p.add_argument("--out", required=True)
            p.add_argument("--plot")
            p.add_argument("--g-max", type=float, default=18.5)
            p.add_argument("--max-sources", type=int, default=800)
            p.add_argument("--aperture-coeff", type=float)
            p.add_argument("--workers", type=int, default=16, help="processes to split the sources across")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    if args.command == "calibrate":
        calibrate(args.seq_dir)
    elif args.command == "solve":
        solve(args.seq_dir)
    else:
        chip = args.chip
        if chip is None:
            primary = fits.getheader(science_files(args.seq_dir)[0], 0)
            planet = planet_parameters(args.planet, float(primary["MJD-OBS"]) + 2400000.5)
            found = chip_containing(args.seq_dir, SkyCoord(planet["ra"] * u.deg, planet["dec"] * u.deg))
            if found is None:
                raise SystemExit(f"{args.planet} is not on any solved chip")
            chip = found[0]
            log.info("%s is on chip %d at (%.0f, %.0f)", args.planet, *found)
        extract(args.seq_dir, chip, args.planet, args.out, plot_path=args.plot, g_max=args.g_max,
                max_sources=args.max_sources, aperture_coeff=args.aperture_coeff, workers=args.workers)


if __name__ == "__main__":
    main()
