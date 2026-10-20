"""Plate solve a frame by matching its stars to Gaia DR3.

Raw archive frames (ESO WFI, MicroObservatory) carry no celestial WCS, and
uploading to nova.astrometry.net makes the image public, so this solves
locally:

1. detect the brightest stars in the frame,
2. project the brightest Gaia stars around the pointing onto a tangent
   plane at roughly the frame's pixel scale,
3. let astroalign match the two point patterns (it finds rotation, scale and
   shift, but not a mirror flip, so both parities are tried),
4. fit a WCS to the matched pairs, then re-match every detection within a
   few pixels and refit.

The Gaia table is returned too, so the solved stars can be photometered and
keyed by Gaia ``source_id``.
"""

import logging

import astroalign as aa
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.wcs.utils import fit_wcs_from_points
from photutils.detection import DAOStarFinder
from scipy.spatial import cKDTree

log = logging.getLogger(__name__)


def gaia_stars(ra, dec, radius_deg, g_max):
    """Gaia DR3 stars within ``radius_deg`` brighter than ``g_max``, brightest first."""
    from astroquery.gaia import Gaia

    query = f"""SELECT source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))
        AND phot_g_mean_mag < {g_max}
        ORDER BY phot_g_mean_mag"""
    return Gaia.launch_job_async(query).get_results()


def detect_stars(data, fwhm=3.0, n_max=200, threshold_sigma=5.0):
    """Pixel positions of the ``n_max`` brightest point sources, brightest first."""
    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    found = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)(data - median)
    if found is None:
        return np.empty((0, 2))
    found.sort("flux", reverse=True)
    found = found[:n_max]
    return np.column_stack([found["xcentroid"], found["ycentroid"]])


def _tangent_plane(ra, dec, ra0, dec0):
    """Gnomonic projection of (ra, dec) about (ra0, dec0), in degrees."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cos_c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cos_c
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cos_c
    return np.degrees(xi), np.degrees(eta)


def _match(xy, sky_xy, max_px):
    """Index pairs (i in xy, j in sky_xy) of mutual nearest neighbours within ``max_px``."""
    tree = cKDTree(sky_xy)
    dist, j = tree.query(xy, distance_upper_bound=max_px)
    ok = np.isfinite(dist)
    i = np.nonzero(ok)[0]
    j = j[ok]
    # Drop Gaia stars claimed by more than one detection.
    unique, counts = np.unique(j, return_counts=True)
    keep = np.isin(j, unique[counts == 1])
    return i[keep], j[keep]


def solve_wcs(data, ra0, dec0, pixel_scale_arcsec, *, search_radius_deg, g_max=15.0,
              n_detect=60, n_gaia=120, fwhm=3.0, sip_degree=2, match_px=3.0, catalogue=None):
    """Fit a celestial WCS to a frame.

    ``ra0, dec0`` is the approximate field centre (e.g. the header pointing)
    and ``pixel_scale_arcsec`` a rough plate scale; astroalign fits the
    exact scale and rotation. Returns ``(wcs, gaia, stats)``: the WCS, the
    Gaia table it was solved against, and a dict with the number of matched
    stars and the RMS residual in arcsec.

    ``catalogue`` is an already-fetched Gaia table (brightest first) covering
    the search circle; when given, no Gaia query is made. Pass one when
    solving many chips or trying many centres over the same field.
    """
    if catalogue is None:
        gaia = gaia_stars(ra0, dec0, search_radius_deg, g_max)
    else:
        sep = SkyCoord(np.asarray(catalogue["ra"]) * u.deg, np.asarray(catalogue["dec"]) * u.deg).separation(
            SkyCoord(ra0 * u.deg, dec0 * u.deg))
        gaia = catalogue[(sep.deg < search_radius_deg) & (np.asarray(catalogue["phot_g_mean_mag"]) < g_max)]
    detected = detect_stars(data, fwhm=fwhm, n_max=max(n_detect, 200))
    if len(detected) < 10 or len(gaia) < 10:
        raise RuntimeError(f"too few stars to solve: {len(detected)} detected, {len(gaia)} Gaia")

    xi, eta = _tangent_plane(np.asarray(gaia["ra"]), np.asarray(gaia["dec"]), ra0, dec0)
    plane = np.column_stack([xi, eta]) * 3600 / pixel_scale_arcsec  # roughly pixels

    best = None
    for parity in (1, -1):
        reference = plane[:n_gaia] * [parity, 1]
        try:
            transform, (src, dst) = aa.find_transform(detected[:n_detect], reference)
        except aa.MaxIterError:
            continue
        if best is None or len(src) > best[2]:
            best = (parity, transform, len(src))
    if best is None:
        raise RuntimeError("astroalign could not match the frame to Gaia")
    parity, transform, n_initial = best
    log.info("pattern match: parity %+d, %d stars, scale %.3f, rotation %.1f deg",
             parity, n_initial, transform.scale, np.degrees(transform.rotation))

    gaia_sky = SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg)

    # First WCS from every detection the similarity transform puts on a Gaia star.
    projected = transform(detected)
    i, j = _match(projected, plane * [parity, 1], match_px)
    wcs = fit_wcs_from_points((detected[i, 0], detected[i, 1]), gaia_sky[j], projection="TAN")

    # Refine: re-match all detections through the fitted WCS, then refit with SIP.
    for _ in range(2):
        gx, gy = wcs.world_to_pixel(gaia_sky)
        on = np.isfinite(gx) & np.isfinite(gy)
        i, j_on = _match(detected, np.column_stack([gx[on], gy[on]]), match_px)
        j = np.nonzero(on)[0][j_on]
        wcs = fit_wcs_from_points((detected[i, 0], detected[i, 1]), gaia_sky[j], projection="TAN",
                                  sip_degree=sip_degree if len(i) >= 30 else None)

    fitted = wcs.pixel_to_world(detected[i, 0], detected[i, 1])
    rms = float(np.sqrt(np.mean(fitted.separation(gaia_sky[j]).arcsec ** 2)))
    stats = {"n_matched": int(len(i)), "rms_arcsec": rms, "parity": parity,
             "pixel_scale_arcsec": float(pixel_scale_arcsec / transform.scale)}
    log.info("WCS: %d stars matched, rms %.2f arcsec, scale %.3f arcsec/px",
             stats["n_matched"], rms, stats["pixel_scale_arcsec"])
    return wcs, gaia, stats


# Search radius, detections and Gaia stars to try, in order. Pattern matching
# needs the two lists to hold roughly the same stars, so these trade a tight
# Gaia circle (header pointing good) against a wide one (pointing poor).
SOLVE_ATTEMPTS = ((0.55, 40, 60), (0.55, 60, 90), (0.45, 50, 60), (0.7, 80, 160))


def solve_wcs_robust(data, centres, pixel_scale_arcsec, *, g_max=15.0, min_matched=20, max_rms_arcsec=2.0,
                     attempts=SOLVE_ATTEMPTS, fwhm=3.0, match_px=3.0, catalogue=None):
    """``solve_wcs`` over several guessed centres and settings; first acceptable solution wins.

    ``centres`` is a list of (ra, dec) guesses, e.g. the header pointing and
    the target's catalogue position. A solution is acceptable with at least
    ``min_matched`` stars and an RMS under ``max_rms_arcsec``.
    """
    errors = []
    for ra0, dec0 in centres:
        for radius, n_detect, n_gaia in attempts:
            try:
                wcs, gaia, stats = solve_wcs(data, ra0, dec0, pixel_scale_arcsec, search_radius_deg=radius,
                                             g_max=g_max, n_detect=n_detect, n_gaia=n_gaia, fwhm=fwhm,
                                             match_px=match_px, catalogue=catalogue)
            except Exception as exc:  # astroalign and fitting failures alike
                errors.append(f"({ra0:.3f}, {dec0:.3f}, r={radius}): {exc}")
                continue
            if stats["n_matched"] >= min_matched and stats["rms_arcsec"] <= max_rms_arcsec:
                stats.update(centre_guess=[float(ra0), float(dec0)], search_radius_deg=radius)
                return wcs, gaia, stats
            errors.append(f"({ra0:.3f}, {dec0:.3f}, r={radius}): poor fit {stats}")
    raise RuntimeError("no acceptable plate solution:\n  " + "\n  ".join(errors))
