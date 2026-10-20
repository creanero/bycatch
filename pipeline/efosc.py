"""ESO 3.6 m NTT EFOSC2 (imaging): calibrate and extract one night.

Usage
-----
    python -m pipeline.efosc SEQ_DIR --planet "WASP-4 b" --out OUT.json [--plot OUT.png] [--aperture-coeff C]

``SEQ_DIR`` is a directory written by ``python -m access.fetch``.

EFOSC2 transit sequences are often read out as a window of the single
2048 x 2048 LORAL CCD (e.g. 700 x 900 px at 0.12"/px), with no pre- or
overscan, while the night's biases and flats are full frames (2060 x 2060:
the imaging area plus 12 prescan columns on the right and 12 overscan rows
on the top). The master bias and flat are built from the full frames and
the science window is cut from them at ``ESO DET WIN1 STRX/STRY``. The
header WCS refers to the full frame, so its reference pixel is moved by
the window offset and the remaining shift is measured against Gaia.

With a field of ~1.5', there may be only one usable comparison star.
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from exotic.extract.apertures import FwhmScaled
from exotic.extract.sequence import Observatory, extract_sequence

from . import store
from .astrometry import detect_stars, gaia_stars
from .lightcurve import LIMB_DARKENING, choose_sources, fit_with_airmass, pick_comparisons
from .wfi import _for_exotic, planet_parameters

log = logging.getLogger("pipeline.efosc")

IMAGING = 2048          # imaging pixels per axis in a full frame
SATURATION = 55000.0    # ADU after bias subtraction


def _manifest(seq_dir):
    return json.loads((Path(seq_dir) / "manifest.json").read_text())


def science_files(seq_dir):
    files = [Path(seq_dir) / "science" / f"{s['dp_id']}.fits" for s in _manifest(seq_dir)["science"]]
    return [p for p in files if p.exists()]


def _full_frames(seq_dir, dp_type):
    out = []
    for c in _manifest(seq_dir)["calibrations"]:
        path = Path(seq_dir) / "calib" / f"{c['dp_id']}.fits"
        if c["dp_type"] == dp_type and path.exists():
            header = fits.getheader(path)
            if header["NAXIS1"] == IMAGING + 12 and header.get("HIERARCH ESO DET WIN1 BINX") == 1:
                out.append(fits.getdata(path).astype(float)[:IMAGING, :IMAGING])
    return out


def masters(seq_dir):
    """Full-frame master bias and normalised flat (imaging area only)."""
    bias = np.median(_full_frames(seq_dir, "BIAS"), axis=0)
    flats = []
    for frame in _full_frames(seq_dir, "SKY,FLAT"):
        level = np.median(frame - bias)
        if 3000 < level < 45000:
            flats.append((frame - bias) / level)
    if len(flats) < 3:
        raise RuntimeError(f"only {len(flats)} usable sky flats")
    flat = np.median(flats, axis=0)
    log.info("master bias from %d frames (level %.1f ADU), flat from %d sky flats",
             len(_full_frames(seq_dir, "BIAS")), np.median(bias), len(flats))
    return bias, flat / np.median(flat), len(flats)


def window_of(header):
    """(row, column) slices of the full-frame imaging area this science window covers."""
    x0 = int(header["HIERARCH ESO DET WIN1 STRX"]) - 1
    y0 = int(header["HIERARCH ESO DET WIN1 STRY"]) - 1
    return slice(y0, y0 + int(header["NAXIS2"])), slice(x0, x0 + int(header["NAXIS1"]))


def calibrated(path, bias, flat):
    data = fits.getdata(path).astype(float)
    header = fits.getheader(path)
    rows, cols = window_of(header)
    return (data - bias[rows, cols]) / flat[rows, cols], header


def window_wcs(header):
    """The header's full-frame WCS with its reference pixel moved to window coordinates."""
    header = header.copy()
    header["CRPIX1"] -= int(header["HIERARCH ESO DET WIN1 STRX"]) - 1
    header["CRPIX2"] -= int(header["HIERARCH ESO DET WIN1 STRY"]) - 1
    return WCS(header)


def refine_shift(image, wcs, gaia, fwhm, tol=4.0):
    """Correct the WCS for a pixel shift by matching detections to Gaia.

    Every detection-Gaia pair votes for an offset; the offset most pairs
    agree with (within ``tol`` px) wins. Works with only a handful of stars.
    """
    det = detect_stars(image, fwhm=fwhm, n_max=30, threshold_sigma=10)
    gx, gy = wcs.world_to_pixel_values(np.asarray(gaia["ra"]), np.asarray(gaia["dec"]))
    offsets = (det[:, None, :] - np.stack([gx, gy], axis=1)[None, :, :]).reshape(-1, 2)
    votes = [np.sum(np.hypot(*(offsets - o).T) < tol) for o in offsets]
    best = offsets[int(np.argmax(votes))]
    agree = offsets[np.hypot(*(offsets - best).T) < tol]
    shift = agree.mean(axis=0)
    if len(agree) < 2:
        raise RuntimeError("could not match the frame to Gaia")
    wcs = wcs.deepcopy()
    wcs.wcs.crpix = wcs.wcs.crpix + shift
    log.info("WCS shift (%.1f, %.1f) px from %d matched stars", shift[0], shift[1], len(agree))
    return wcs, {"shift_px": shift.tolist(), "n_matched": int(len(agree))}


def profile_fwhm(image, x, y, max_r=60):
    """FWHM (px) from a star's azimuthally averaged radial profile.

    EFOSC2 transit sequences are usually defocused into broad, flat-topped
    images, which EXOTIC's Gaussian FWHM measurement can't fit. Twice the
    radius at which the sky-subtracted profile falls to half its central
    value is robust to the shape.
    """
    yv, xv = np.mgrid[:image.shape[0], :image.shape[1]]
    r = np.hypot(xv - x, yv - y)
    sky = np.median(image[(r > max_r) & (r < max_r + 20)])
    radii = np.arange(0, max_r)
    profile = np.array([np.median(image[(r >= a) & (r < a + 1)]) for a in radii]) - sky
    half = profile[:3].mean() / 2
    below = np.nonzero(profile < half)[0]
    return 2.0 * float(radii[below[0]]) if len(below) else float("nan")


def track_offsets(files, bias, flat, fwhm, start, max_step=50.0):
    """Per-frame (dx, dy) of the star that starts at ``start`` = (x, y), relative to the first frame.

    Follows the detection nearest the star's last position (not the
    brightest: in a small field that can swap between the target and a
    comparison of similar brightness). Frames with nothing found nearby
    reuse the previous offset.
    """
    start = np.asarray(start, dtype=float)
    position = start.copy()
    offsets = []
    for path in files:
        image, _ = calibrated(path, bias, flat)
        found = detect_stars(np.nan_to_num(image), fwhm=fwhm, n_max=10, threshold_sigma=10)
        if len(found):
            nearest = found[np.argmin(np.hypot(*(found - position).T))]
            if np.hypot(*(nearest - position)) < max_step:
                position = nearest
        offsets.append(position - start)
    return np.array(offsets)


def shifted(image, offset):
    """``image`` moved back by the whole-pixel part of ``offset``; uncovered edges set to the median sky.

    Whole-pixel shifts move pixels without resampling, so photometry is
    unchanged; the sub-pixel remainder is left to EXOTIC's centroid fit.
    """
    dx, dy = (int(round(v)) for v in offset)
    out = np.full_like(image, np.nanmedian(image))
    ny, nx = image.shape
    out[max(0, -dy):ny - max(0, dy), max(0, -dx):nx - max(0, dx)] = \
        image[max(0, dy):ny - max(0, -dy), max(0, dx):nx - max(0, -dx)]
    return out


def run(seq_dir, planet_name, out_path, plot_path=None, aperture_coeff=None, g_max=19.0, fwhm=20.0):
    seq_dir = Path(seq_dir)
    files = science_files(seq_dir)
    bias, flat, n_flats = masters(seq_dir)
    reference, header = calibrated(files[0], bias, flat)
    reference = np.nan_to_num(reference)
    obs = Observatory(latitude=float(header["HIERARCH ESO TEL GEOLAT"]),
                      longitude=float(header["HIERARCH ESO TEL GEOLON"]),
                      elevation=float(header["HIERARCH ESO TEL GEOELEV"]))
    planet = planet_parameters(planet_name, float(header["MJD-OBS"]) + 2400000.5)
    target = SkyCoord(planet["ra"] * u.deg, planet["dec"] * u.deg)

    wcs = window_wcs(header)
    centre = wcs.pixel_to_world(reference.shape[1] / 2, reference.shape[0] / 2)
    gaia = gaia_stars(centre.ra.deg, centre.dec.deg, 0.05, g_max + 2.5)
    wcs, astrometry = refine_shift(reference, wcs, gaia, fwhm)
    picked = choose_sources(reference, wcs, gaia, target, g_max, edge_px=40, blend_px=40,
                            saturation=SATURATION)
    gaia = gaia[picked]
    log.info("photometering %d sources: G = %s", len(gaia), np.round(np.asarray(gaia["phot_g_mean_mag"]), 1).tolist())

    # Too few stars for astroalign: register by tracking the target and shifting whole pixels.
    tx, ty = wcs.world_to_pixel(target)
    offsets = track_offsets(files, bias, flat, fwhm, (float(tx), float(ty)))
    drift = float(np.max(np.hypot(*offsets.T)))
    log.info("target drifts up to %.1f px; frames shifted back by whole pixels", drift)

    def frames():
        for i, path in enumerate(files):
            image, h = calibrated(path, bias, flat)
            if i == 0:
                h = h.copy()
                for key in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                            "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"):
                    h.remove(key, ignore_missing=True)
                h.update(wcs.to_header())
            else:
                h = h.copy()
                for key in ("CTYPE1", "CTYPE2"):  # stale full-frame WCS; RA/Dec come via the reference
                    h.remove(key, ignore_missing=True)
            yield _for_exotic(shifted(image, offsets[i])), h

    gain = float(header["HIERARCH ESO DET OUT1 CONAD"])
    read_noise = float(header["HIERARCH ESO DET OUT1 RON"])
    measured_fwhm = profile_fwhm(reference, float(tx), float(ty))
    scaled = FwhmScaled() if aperture_coeff is None else FwhmScaled(aperture_coeff=aperture_coeff)
    apertures = scaled.apertures_for(measured_fwhm)
    log.info("FWHM %.1f px from the target's radial profile; aperture r=%.1f, annulus %.1f-%.1f px",
             measured_fwhm, apertures.aperture_r, apertures.annulus_r_in, apertures.annulus_r_out)
    table = extract_sequence(frames(), sky_positions=list(zip(np.asarray(gaia["ra"], float),
                                                              np.asarray(gaia["dec"], float))),
                             ra=target.ra.deg, dec=target.dec.deg, observatory=obs, apertures=apertures,
                             gain=gain, read_noise=read_noise, align=False)
    grid = table.reshape(len(files), len(gaia))
    flux = grid["flux"]
    bjd, airmass = grid["bjd_tdb"][:, 0], grid["airmass"][:, 0]
    good = grid["aligned"][:, 0] & np.isfinite(flux[:, 0]) & (flux[:, 0] > 0)

    comps, comp_scatter = pick_comparisons(flux[good][:, 1:], flux[good][:, 0], n_max=5, ratio_range=(0.2, 5.0))
    comps = [c + 1 for c in comps]
    relative = flux[:, 0] / flux[:, comps].sum(axis=1)
    relative = relative / np.nanmedian(relative[good])
    values, errors, extra = fit_with_airmass(bjd[good], relative[good], airmass[good], planet)
    fit = {"tmid_bjd": values["tmid"], "tmid_err_min": errors["tmid"] * 1440, "rprs": values["rprs"],
           "rprs_err": errors["rprs"], "scale": values["scale"],
           "o_minus_c_min": (values["tmid"] - planet["tmid_predicted"]) * 1440,
           "airmass_poly": extra["airmass_poly"], "residual_rms": extra["rms"],
           "n_points": int(extra["kept"].sum()), "n_clipped": extra["n_clipped"]}
    log.info("fit: T_mid %.5f +- %.2f min, O-C %+.2f min (prediction +- %.2f), Rp/Rs %.4f +- %.4f, rms %.4f",
             fit["tmid_bjd"], fit["tmid_err_min"], fit["o_minus_c_min"], planet["tmid_predicted_sigma_min"],
             fit["rprs"], fit["rprs_err"], fit["residual_rms"])

    roles = ["target"] + ["comparison" if j in comps else "field" for j in range(1, len(gaia))]
    sources = [{"gaia_dr3": int(gaia["source_id"][j]), "g_mag": float(gaia["phot_g_mean_mag"][j]),
                "ra": float(gaia["ra"][j]), "dec": float(gaia["dec"][j]), "role": roles[j],
                "flux": grid["flux"][:, j], "uncertainty": grid["uncertainty"][:, j],
                "x": grid["x_centroid"][:, j], "y": grid["y_centroid"][:, j]} for j in range(len(gaia))]
    manifest = _manifest(seq_dir)
    sequence = {"archive": "ESO", "instrument": "EFOSC", "object": header.get("OBJECT"),
                "filter": header.get("HIERARCH ESO INS FILT1 NAME"), "exposure_s": float(header["EXPTIME"]),
                "night": manifest["target"].get("night"), "dp_ids": [p.stem for p in files],
                "window": {"strx": int(header["HIERARCH ESO DET WIN1 STRX"]),
                           "stry": int(header["HIERARCH ESO DET WIN1 STRY"]),
                           "shape": list(reference.shape)},
                "site": {"latitude": obs.latitude, "longitude": obs.longitude, "elevation": obs.elevation},
                "apertures": {k: float(grid[k][0, 0]) for k in ("aperture_r", "annulus_r_in", "annulus_r_out")},
                "gain": gain, "read_noise": read_noise, "uncertainties": "header gain/read noise",
                "calibration": {"bias_level_adu": float(np.median(bias)), "n_flats": n_flats,
                                "flat_window_alignment": "nominal STRX/STRY offset; not independently verified"},
                "astrometry": astrometry, "registration": {"method": "target tracking, whole-pixel shifts", "max_drift_px": drift,
                                 "offsets_px": offsets.tolist()},
                "fwhm_px": measured_fwhm, "g_max": g_max, "limb_darkening": LIMB_DARKENING}
    frames_out = {"bjd_tdb": bjd, "jd_utc": grid["jd_utc"][:, 0], "airmass": airmass,
                  "aligned": grid["aligned"][:, 0], "used": good}
    analysis = {"planet": planet, "comparisons": {"gaia_dr3": [int(gaia["source_id"][c]) for c in comps],
                                                  "point_to_point_scatter": list(comp_scatter.values())},
                "relative_flux": relative, "transit_fit": fit}
    store.write(out_path, sequence, frames_out, sources, analysis)
    log.info("wrote %s", out_path)
    if plot_path:
        from .refit import refit
        refit(out_path, plot_path)  # same fit; refit draws the standard plot
    return fit


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_dir")
    parser.add_argument("--planet", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--plot")
    parser.add_argument("--aperture-coeff", type=float)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    run(args.seq_dir, args.planet, args.out, plot_path=args.plot, aperture_coeff=args.aperture_coeff)


if __name__ == "__main__":
    main()
