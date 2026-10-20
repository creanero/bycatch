"""Known transiting planets, and which archive sequences catch a transit.

Ephemerides come from the NASA Exoplanet Archive's ``pscomppars`` table.
Transit times there are BJD_TDB; sequence times are MJD (UTC). The few
minutes of barycentric and UTC/TDB difference don't matter for deciding
whether a transit falls in a night, so they are ignored here.
"""

import numpy as np
import pyvo
from astropy.table import Table

NEA_TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP"

# Used when the archive has no transit duration for a planet.
DEFAULT_DURATION_HOURS = 3.0

# Half-width of each instrument's field of view, in degrees.
FIELD_HALF_WIDTH = {"WFI": 0.27, "OMEGACAM": 0.48, "FORS1": 0.057, "FORS2": 0.057, "EFOSC": 0.034,
                    "HAWKI": 0.06, "SOFI": 0.041}


def _value(v, default=0.0):
    """A catalogue value as a float, with ``default`` for missing (masked) or non-finite entries."""
    if np.ma.is_masked(v):
        return default
    v = float(v)
    return v if np.isfinite(v) else default


def transiting_planets(max_dec=30.0):
    """Transiting planets with a period and mid-transit time."""
    query = f"""SELECT pl_name, ra, dec, pl_orbper, pl_orbpererr1, pl_tranmid, pl_tranmiderr1,
        pl_trandur, pl_trandep, pl_ratror, sy_vmag
        FROM pscomppars WHERE tran_flag=1 AND pl_tranmid IS NOT NULL AND pl_orbper IS NOT NULL
        AND dec < {max_dec}"""
    return pyvo.dal.TAPService(NEA_TAP_URL).search(query).to_table()


def transits_in(sequences, planets, half_width=FIELD_HALF_WIDTH):
    """Transits of known planets that fall inside a sequence's field and time window.

    ``sequences`` is the table from ``eso.find_sequences``. Returns one row per
    (sequence, transit), with ``coverage`` the fraction of the transit
    observed, ``before_h`` and ``after_h`` the hours of baseline either side
    (negative if ingress or egress is missed), and ``tc_sigma_min`` the
    propagated 1-sigma uncertainty of the predicted mid-transit time.
    """
    rows = []
    for seq in sequences:
        h = half_width[seq["instrument"]]
        dra = ((planets["ra"] - seq["ra"]) * np.cos(np.radians(seq["dec"])) + 180) % 360 - 180
        in_field = (np.abs(dra) < h) & (np.abs(planets["dec"] - seq["dec"]) < h)
        start, end = seq["mjd_start"] + 2400000.5, seq["mjd_end"] + 2400000.5
        for pl in planets[in_field]:
            period, tc0 = pl["pl_orbper"], pl["pl_tranmid"]
            hours = _value(pl["pl_trandur"], DEFAULT_DURATION_HOURS)
            duration = hours / 24
            for n in np.arange(np.ceil((start - tc0) / period) - 1, np.floor((end - tc0) / period) + 2):
                tc = tc0 + n * period
                ingress, egress = tc - duration / 2, tc + duration / 2
                coverage = (min(egress, end) - max(ingress, start)) / duration
                if coverage <= 0:
                    continue
                sigma = np.hypot(_value(pl["pl_tranmiderr1"]), abs(n) * _value(pl["pl_orbpererr1"])) * 24 * 60
                rows.append((pl["pl_name"], seq["instrument"], seq["prog_id"], seq["object"], seq["filter"],
                             int(seq["night"]), int(seq["n_frames"]), tc, min(coverage, 1.0),
                             (ingress - start) * 24, (end - egress) * 24, sigma,
                             pl["pl_trandep"], pl["sy_vmag"]))
    names = ("planet", "instrument", "prog_id", "object", "filter", "night", "n_frames", "tc_bjd",
             "coverage", "before_h", "after_h", "tc_sigma_min", "depth_pct", "vmag")
    table = Table(rows=rows, names=names) if rows else Table(names=names)
    if len(table):
        table.sort("coverage", reverse=True)
    return table
