"""Extract and fit a MicroObservatory night with EXOTIC's extraction API.

Usage
-----
    python -m pipeline.microobs FRAMES_DIR --inits INITS.json --out RESULTS.json [--plot LC.png]

Steps:

1. Read the frames (darks, named ``Dark*``, are applied if present).
2. Plate solve the first frame against Gaia DR3 (``pipeline.astrometry``);
   MicroObservatory headers have no WCS and the header pointing can be
   ~0.3 deg off.
3. Choose sources: the target, plus every unblended, unsaturated Gaia star
   on the frame brighter than ``--g-max``.
4. ``exotic.extract.sequence.extract_sequence`` over all frames, with the
   sources given as RA/Dec through the solved WCS.
5. Pick comparison stars from the field (similar brightness, lowest
   scatter against each other), form the target's relative light curve,
   and fit it with ``exotic.api.quick_fit.fit_transit``, with a linear
   airmass trend divided out between iterations.
6. Write every source's light curve, the fit and the run metadata to JSON
   (``pipeline.store``).

Gain and read noise for MicroObservatory are not calibrated yet, so the
per-point uncertainties are placeholders and the fit scales its errors from
the residuals.
"""

import argparse
import json
import logging
from datetime import timezone
from pathlib import Path

import dateutil.parser
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time

from exotic.api.elca import transit
from exotic.api.quick_fit import fit_transit
from exotic.extract.apertures import FwhmScaled
from exotic.extract.sequence import Observatory, extract_sequence

from . import store
from .astrometry import gaia_stars, solve_wcs_robust

log = logging.getLogger("pipeline.microobs")

PIXEL_SCALE_GUESS = 5.2       # arcsec/px; the solve fits the real value
GAIN, READ_NOISE = 1.0, 0.0   # placeholders, not a calibration
# Rough Claret 4-parameter coefficients for a Sun-like star in a clear
# filter, as in EXOTIC's examples. Good enough for timing, not for depth.
LIMB_DARKENING = (0.0, 0.9, 0.0, -0.2)


def read_frames(frames_dir):
    """Science and dark file lists, science frames in time order."""
    files = sorted(p for p in Path(frames_dir).iterdir() if p.suffix.lower() in (".fits", ".fit", ".fts"))
    darks = [p for p in files if p.name.lower().startswith("dark")]
    science = [p for p in files if p not in darks]
    science.sort(key=lambda p: float(fits.getheader(p)["MJD-OBS"]))
    return science, darks


def master_dark(darks):
    if not darks:
        return None
    return np.median([fits.getdata(p).astype(float) for p in darks], axis=0)


def date_obs_jd_mid(header):
    """JD (UTC) of mid-exposure from DATE-OBS, honouring its timezone offset."""
    when = dateutil.parser.parse(header["DATE-OBS"]).astimezone(timezone.utc).replace(tzinfo=None)
    return Time(when, scale="utc").jd + float(header.get("EXPTIME", 0.0)) / 2 / 86400


def choose_sources(reference, wcs, gaia, target, g_max, edge_px=20, blend_px=12, blend_dmag=2.5):
    """Indices into ``gaia`` of usable field stars, and the target's index.

    Drops stars within ``edge_px`` of the frame edge, stars with a neighbour
    no more than ``blend_dmag`` fainter within ``blend_px``, and stars whose
    peak reaches 95% of the frame's maximum (saturated).
    """
    sky = SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg)
    target_index = int(np.argmin(sky.separation(target)))
    if sky[target_index].separation(target) > 5 * u.arcsec:
        raise RuntimeError("target has no Gaia counterpart within 5 arcsec")

    x, y = wcs.world_to_pixel(sky)
    ny, nx = reference.shape
    g = np.asarray(gaia["phot_g_mean_mag"], dtype=float)
    on = (x > edge_px) & (x < nx - edge_px) & (y > edge_px) & (y < ny - edge_px)
    saturation = 0.95 * reference.max()

    keep = []
    for k in np.nonzero(on & (g < g_max))[0]:
        d = np.hypot(x - x[k], y - y[k])
        if np.any((d > 0) & (d < blend_px) & (g < g[k] + blend_dmag)):
            continue
        xi, yi = int(round(x[k])), int(round(y[k]))
        if reference[max(yi - 3, 0):yi + 4, max(xi - 3, 0):xi + 4].max() >= saturation:
            continue
        keep.append(int(k))
    keep = [k for k in keep if k != target_index]
    return [target_index] + keep


def pick_comparisons(flux, target_flux, n_max=8, ratio_range=(0.3, 5.0)):
    """Comparison-star columns of ``flux`` (frames x candidates).

    Candidates must be finite in every frame and within ``ratio_range`` of
    the target's median flux. The ensemble is then pruned, dropping the star
    that scatters most against the rest, until ``n_max`` remain.
    """
    median_t = np.nanmedian(target_flux)
    med = np.nanmedian(flux, axis=0)
    ok = np.all(np.isfinite(flux), axis=0) & (med > ratio_range[0] * median_t) & (med < ratio_range[1] * median_t)
    chosen = list(np.nonzero(ok)[0])
    if len(chosen) < 2:
        raise RuntimeError("fewer than two usable comparison stars")

    def scatter(j, others):
        ratio = flux[:, j] / flux[:, others].sum(axis=1)
        ratio = ratio / np.median(ratio)
        return np.std(np.diff(ratio)) / np.sqrt(2)  # point-to-point, insensitive to slow trends

    while len(chosen) > n_max:
        worst = max(chosen, key=lambda j: scatter(j, [c for c in chosen if c != j]))
        chosen.remove(worst)
    return chosen, {int(j): float(scatter(j, [c for c in chosen if c != j])) for j in chosen}


def fit_with_airmass(times, flux, airmass, planet, n_iter=3):
    """``fit_transit`` with a linear airmass trend divided out between fits."""
    tmid_guess = planet["tmid_predicted"]
    trend = np.ones_like(flux)
    for _ in range(n_iter):
        values, errors = fit_transit(times, flux / trend, planet["per"], planet["ars"], planet["inc"],
                                     tmid_guess, planet["rprs"], LIMB_DARKENING,
                                     ecc=planet["ecc"], omega=planet["omega"])
        model = values["scale"] * transit(times, dict(
            per=planet["per"], ars=planet["ars"], inc=planet["inc"], ecc=planet["ecc"], omega=planet["omega"],
            u0=LIMB_DARKENING[0], u1=LIMB_DARKENING[1], u2=LIMB_DARKENING[2], u3=LIMB_DARKENING[3],
            tmid=values["tmid"], rprs=values["rprs"]))
        coeffs = np.polyfit(airmass - 1.0, flux / model, 1)
        trend = np.polyval(coeffs, airmass - 1.0)
        tmid_guess = values["tmid"]
    residuals = flux / trend - model
    return values, errors, {"airmass_poly": coeffs.tolist(), "rms": float(np.std(residuals)),
                            "model": model, "detrended": flux / trend}


def planet_from_inits(inits, bjd_mid):
    pp = inits["planetary_parameters"]
    per = float(pp["Orbital Period (days)"])
    t0 = float(pp["Published Mid-Transit Time (BJD-UTC)"])
    epoch = int(np.round((bjd_mid - t0) / per))
    sigma = np.hypot(float(pp["Mid-Transit Time Uncertainty"]), abs(epoch) * float(pp["Orbital Period Uncertainty"]))
    return {"name": pp["Planet Name"], "ra": float(pp["Target Star RA"]), "dec": float(pp["Target Star Dec"]),
            "per": per, "ars": float(pp["Ratio of Distance to Stellar Radius (a/Rs)"]),
            "inc": float(pp["Orbital Inclination (deg)"]), "ecc": float(pp["Orbital Eccentricity (0 if null)"]),
            "omega": float(pp["Argument of Periastron (deg)"]),
            "rprs": float(pp["Ratio of Planet to Stellar Radius (Rp/Rs)"]),
            "t0_published": t0, "epoch": epoch, "tmid_predicted": t0 + epoch * per,
            "tmid_predicted_sigma_min": float(sigma * 1440)}


def run(frames_dir, inits_path, out_path, g_max=13.5, plot_path=None, aperture_coeff=None):
    inits = json.loads(Path(inits_path).read_text())
    ui = inits["user_info"]
    site = Observatory(latitude=float(ui["Obs. Latitude"]), longitude=float(ui["Obs. Longitude"]),
                       elevation=float(ui["Obs. Elevation (meters; Note: leave blank if unknown)"]))
    pp = inits["planetary_parameters"]
    target = SkyCoord(float(pp["Target Star RA"]) * u.deg, float(pp["Target Star Dec"]) * u.deg)

    science, darks = read_frames(frames_dir)
    dark = master_dark(darks)
    log.info("%d science frames, %d darks", len(science), len(darks))

    ref_header = fits.getheader(science[0])
    reference = fits.getdata(science[0]).astype(float)
    if dark is not None:
        reference = reference - dark
    wcs, _, solve_stats = solve_wcs_robust(
        reference, [(ref_header["RA"], ref_header["DEC"]), (target.ra.deg, target.dec.deg)], PIXEL_SCALE_GUESS)

    # Gaia over the whole solved frame, not just the circle the solve used.
    ny, nx = reference.shape
    centre = wcs.pixel_to_world(nx / 2, ny / 2)
    half_diag = np.hypot(nx, ny) / 2 * solve_stats["pixel_scale_arcsec"] / 3600
    gaia = gaia_stars(centre.ra.deg, centre.dec.deg, half_diag * 1.05, g_max + 2.5)
    picked = choose_sources(reference, wcs, gaia, target, g_max)
    gaia = gaia[picked]
    log.info("photometering %d sources (target + %d field stars)", len(gaia), len(gaia) - 1)

    def frames():
        for i, path in enumerate(science):
            header = fits.getheader(path)
            if i == 0:
                header = header.copy()
                header.update(wcs.to_header(relax=True))
            yield fits.getdata(path).astype(float), header

    apertures = FwhmScaled() if aperture_coeff is None else FwhmScaled(aperture_coeff=aperture_coeff)
    sky_positions = list(zip(np.asarray(gaia["ra"], float), np.asarray(gaia["dec"], float)))
    table = extract_sequence(frames(), sky_positions=sky_positions, ra=target.ra.deg, dec=target.dec.deg,
                             observatory=site, apertures=apertures, gain=GAIN, read_noise=READ_NOISE,
                             dark=dark)
    n_frames, n_sources = len(science), len(gaia)
    grid = table.reshape(n_frames, n_sources)
    flux = grid["flux"].copy()
    bjd = grid["bjd_tdb"][:, 0]
    jd_utc = grid["jd_utc"][:, 0]
    airmass = grid["airmass"][:, 0]
    aligned = grid["aligned"][:, 0]

    # EXOTIC times MicroObservatory frames from MJD-OBS, which has 3 decimals
    # (~86 s). DATE-OBS is to the millisecond: measure the difference and use it.
    jd_date_obs = np.array([date_obs_jd_mid(fits.getheader(p)) for p in science])
    timing_offset_s = (jd_date_obs - jd_utc) * 86400
    bjd_precise = bjd + (jd_date_obs - jd_utc)

    good = aligned & np.isfinite(flux[:, 0]) & (flux[:, 0] > 0)
    comps, comp_scatter = pick_comparisons(flux[good][:, 1:], flux[good][:, 0])
    comps = [c + 1 for c in comps]
    ensemble = flux[:, comps].sum(axis=1)
    relative = flux[:, 0] / ensemble
    relative = relative / np.nanmedian(relative[good])

    planet = planet_from_inits(inits, float(np.median(bjd)))
    fits_out = {}
    for label, times in (("mjd_obs", bjd), ("date_obs", bjd_precise)):
        values, errors, extra = fit_with_airmass(times[good], relative[good], airmass[good], planet)
        fits_out[label] = {
            "tmid_bjd": values["tmid"], "tmid_err_min": errors["tmid"] * 1440,
            "rprs": values["rprs"], "rprs_err": errors["rprs"], "scale": values["scale"],
            "o_minus_c_min": (values["tmid"] - planet["tmid_predicted"]) * 1440,
            "airmass_poly": extra["airmass_poly"], "residual_rms": extra["rms"]}
        if label == "date_obs":
            detrended, model = extra["detrended"], extra["model"]

    roles = ["target"] + ["comparison" if j in comps else "field" for j in range(1, n_sources)]
    sources = [{"gaia_dr3": int(gaia["source_id"][j]), "g_mag": float(gaia["phot_g_mean_mag"][j]),
                "ra": float(gaia["ra"][j]), "dec": float(gaia["dec"][j]), "role": roles[j],
                "flux": grid["flux"][:, j], "uncertainty": grid["uncertainty"][:, j],
                "x": grid["x_centroid"][:, j], "y": grid["y_centroid"][:, j]}
               for j in range(n_sources)]
    sequence = {
        "archive": "MicroObservatory", "instrument": ref_header.get("CREATOR"),
        "object": ref_header.get("OBJECT"), "filter": ref_header.get("FILTER"),
        "exposure_s": float(ref_header.get("EXPTIME", 0.0)),
        "files": [p.name for p in science], "darks": [p.name for p in darks],
        "site": {"latitude": site.latitude, "longitude": site.longitude, "elevation": site.elevation},
        "apertures": {k: float(table[k][0]) for k in ("aperture_r", "annulus_r_in", "annulus_r_out")},
        "gain": GAIN, "read_noise": READ_NOISE, "uncertainties": "placeholder gain/read noise",
        "astrometry": solve_stats, "g_max": g_max, "limb_darkening": LIMB_DARKENING,
    }
    frames_out = {"bjd_tdb": bjd, "bjd_tdb_date_obs": bjd_precise, "jd_utc": jd_utc, "airmass": airmass,
                  "aligned": aligned, "used": good}
    analysis = {"planet": planet, "comparisons": {"gaia_dr3": [int(gaia["source_id"][c]) for c in comps],
                                                  "point_to_point_scatter": list(comp_scatter.values())},
                "relative_flux": relative, "transit_fit": fits_out,
                "timing_offset_s": {"median": float(np.median(timing_offset_s)),
                                    "min": float(timing_offset_s.min()), "max": float(timing_offset_s.max())}}
    store.write(out_path, sequence, frames_out, sources, analysis)
    log.info("wrote %s", out_path)

    if plot_path:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = (bjd_precise[good] - planet["tmid_predicted"]) * 24
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(t, detrended, ".", color="0.4", label="WASP-36 / comparison ensemble, airmass-detrended")
        order = np.argsort(t)
        ax.plot(t[order], model[order], "-", color="C3", label="fit_transit")
        f = fits_out["date_obs"]
        ax.axvline(f["o_minus_c_min"] / 60, color="C3", ls=":", lw=1)
        ax.axvline(0, color="0.6", ls="--", lw=1, label="predicted T_mid")
        ax.set_xlabel("hours from predicted mid-transit")
        ax.set_ylabel("relative flux")
        ax.set_title(f"{planet['name']}: T_mid {f['tmid_bjd']:.5f} BJD ± {f['tmid_err_min']:.1f} min, "
                     f"O-C {f['o_minus_c_min']:+.1f} min, Rp/Rs {f['rprs']:.3f}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=120)
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("frames_dir")
    parser.add_argument("--inits", required=True, help="EXOTIC inits JSON with planetary parameters and site")
    parser.add_argument("--out", required=True, help="results JSON path")
    parser.add_argument("--plot", help="light-curve PNG path")
    parser.add_argument("--aperture-coeff", type=float,
                        help="aperture radius in units of the measured FWHM (default: FwhmScaled's)")
    parser.add_argument("--g-max", type=float, default=13.5, help="faintest Gaia G magnitude to photometer")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    run(args.frames_dir, args.inits, args.out, g_max=args.g_max, plot_path=args.plot,
        aperture_coeff=args.aperture_coeff)


if __name__ == "__main__":
    main()
