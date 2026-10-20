"""Instrument-independent light-curve steps: choosing sources and comparison stars, and the transit fit."""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from exotic.api.elca import transit
from exotic.api.quick_fit import fit_transit

# Rough Claret 4-parameter coefficients, as in EXOTIC's examples. Good enough
# for timing, not for depth.
LIMB_DARKENING = (0.0, 0.9, 0.0, -0.2)


def match_target(gaia, target):
    """Index of the target's Gaia counterpart.

    The nearest star within 2 arcsec if there is one (catalogue coordinates).
    Otherwise the brightest within 15 arcsec, because EXOTIC inits
    coordinates are often rounded to 1 s of RA / 1 arcsec, ~10 arcsec at high
    declination.
    """
    sky = SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg)
    sep = sky.separation(target)
    if sep.min() < 2 * u.arcsec:
        return int(np.argmin(sep))
    near = np.nonzero(sep < 15 * u.arcsec)[0]
    if len(near) == 0:
        raise RuntimeError("target has no Gaia counterpart within 15 arcsec")
    return int(near[np.argmin(np.asarray(gaia["phot_g_mean_mag"])[near])])


def choose_sources(reference, wcs, gaia, target, g_max, edge_px=20, blend_px=12, blend_dmag=2.5,
                   saturation=None, max_sources=None):
    """Indices into ``gaia`` of usable field stars, target first.

    Drops stars within ``edge_px`` of the frame edge, stars with a neighbour
    no more than ``blend_dmag`` fainter within ``blend_px``, and stars whose
    peak reaches ``saturation`` (default: 95% of the frame's maximum). With
    ``max_sources``, keeps the brightest that many. ``target`` may be None
    for a frame (chip) without one.
    """
    sky = SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg)
    target_index = None if target is None else match_target(gaia, target)

    x, y = wcs.world_to_pixel(sky)
    ny, nx = reference.shape
    g = np.asarray(gaia["phot_g_mean_mag"], dtype=float)
    on = (x > edge_px) & (x < nx - edge_px) & (y > edge_px) & (y < ny - edge_px)
    if saturation is None:
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
    keep = sorted((k for k in keep if k != target_index), key=lambda k: g[k])
    if max_sources is not None:
        keep = keep[:max_sources]
    return keep if target_index is None else [target_index] + keep


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


def _model(times, planet, tmid, rprs):
    return transit(times, dict(
        per=planet["per"], ars=planet["ars"], inc=planet["inc"], ecc=planet["ecc"], omega=planet["omega"],
        u0=LIMB_DARKENING[0], u1=LIMB_DARKENING[1], u2=LIMB_DARKENING[2], u3=LIMB_DARKENING[3],
        tmid=tmid, rprs=rprs))


def _grid_start(times, flux, planet, tmid0, half_width_days=0.06, n_tmid=121):
    """Best (tmid, rprs) on a coarse grid, with the flux scale solved analytically.

    A starting point for ``fit_transit``. Its least-squares optimiser can sit
    at its initial values when the data are noisy.
    """
    best = (np.inf, tmid0, planet["rprs"])
    for tmid in tmid0 + np.linspace(-half_width_days, half_width_days, n_tmid):
        for rprs in planet["rprs"] * np.array([0.5, 0.75, 1.0, 1.25, 1.5]):
            m = _model(times, planet, tmid, rprs)
            scale = np.dot(flux, m) / np.dot(m, m)
            chi2 = np.sum((flux - scale * m) ** 2)
            if chi2 < best[0]:
                best = (chi2, tmid, rprs)
    return best[1], best[2]


def fit_with_airmass(times, flux, airmass, planet, n_iter=3, clip_sigma=4.0):
    """``fit_transit`` with a linear airmass trend divided out and outliers clipped between fits.

    Times are fitted relative to the integer day before the first point.
    ``fit_transit``'s finite-difference steps scale with the parameter
    value, so at BJD ~2.45e6 the step on ``tmid`` is ~50 minutes and the
    optimiser can fail to move. The start comes from a coarse grid
    (``_grid_start``). Points more than ``clip_sigma`` robust standard
    deviations from the model are dropped from the next iteration;
    ``extra["kept"]`` marks the points used in the final fit.
    """
    times = np.asarray(times, dtype=float)
    flux = np.asarray(flux, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    t_ref = np.floor(times.min())
    t = times - t_ref
    local = dict(planet, tmid_predicted=planet["tmid_predicted"] - t_ref)

    kept = np.isfinite(t) & np.isfinite(flux) & np.isfinite(airmass)
    trend = np.ones_like(flux)
    tmid_guess, rprs_guess = _grid_start(t[kept], flux[kept], local, local["tmid_predicted"])
    for iteration in range(n_iter):
        values, errors = fit_transit(t[kept], (flux / trend)[kept], planet["per"], planet["ars"], planet["inc"],
                                     tmid_guess, rprs_guess, LIMB_DARKENING,
                                     ecc=planet["ecc"], omega=planet["omega"])
        model = values["scale"] * _model(t, local, values["tmid"], values["rprs"])
        coeffs = np.polyfit(airmass[kept] - 1.0, flux[kept] / model[kept], 1)
        trend = np.polyval(coeffs, airmass - 1.0)
        if iteration < n_iter - 1:
            residuals = flux / trend - model
            mad = 1.4826 * np.median(np.abs(residuals[kept] - np.median(residuals[kept])))
            kept = kept & (np.abs(residuals) < clip_sigma * mad)
        tmid_guess, rprs_guess = values["tmid"], abs(values["rprs"])
    values = dict(values, tmid=values["tmid"] + t_ref, rprs=abs(values["rprs"]))
    residuals = flux / trend - model
    return values, errors, {"airmass_poly": coeffs.tolist(), "rms": float(np.std(residuals[kept])),
                            "model": model, "detrended": flux / trend, "kept": kept,
                            "n_clipped": int(np.sum(~kept))}
