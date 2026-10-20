"""Search and download raw frames from the ESO Science Archive.

Metadata comes from the ESO TAP service (table ``dbo.raw``); files come from
the ESO dataportal through ``astroquery.eso``. Public frames need no login.

ESO's ADQL can't group by an expression, so frames are grouped into nights
client side. A night is ``floor(mjd_obs - 0.5)``: it changes at 12:00 UT,
which is morning in Chile, so one observing night is one integer.
"""

import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyvo
from astropy.table import Table

TAP_URL = "https://archive.eso.org/tap_obs"

FRAME_COLUMNS = ("dp_id", "instrument", "prog_id", "object", "dp_cat", "dp_type", "dp_tech",
                 "filter_path", "exposure", "mjd_obs", "ra", "dec", "access_estsize", "release_date")

log = logging.getLogger(__name__)


def night_of(mjd_obs):
    """Integer night label for an MJD (changes at 12:00 UT)."""
    return np.floor(np.asarray(mjd_obs) - 0.5).astype(int)


def _tap():
    return pyvo.dal.TAPService(TAP_URL)


def _quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def query_frames(where, columns=FRAME_COLUMNS, maxrec=1_000_000):
    """Rows of ``dbo.raw`` matching an ADQL ``where`` clause, with a ``night`` column."""
    query = f"SELECT {', '.join(columns)} FROM dbo.raw WHERE {where}"
    table = _tap().run_async(query, maxrec=maxrec).to_table()
    for name in table.colnames:
        if table[name].dtype.kind in "OU":
            table[name] = [str(v) for v in table[name]]
    if "mjd_obs" in table.colnames:
        table["night"] = night_of(table["mjd_obs"])
    return table


def find_sequences(instruments, *, min_frames=40, min_span_days=2 / 24, max_pointing_scatter=0.1,
                   released_before="2026-10-01"):
    """Long single-pointing imaging sequences, one row per (object, filter, night).

    A sequence is the public SCIENCE imaging frames of one instrument,
    programme, object and filter in one night, kept when it has at least
    ``min_frames`` frames spanning at least ``min_span_days``, with the
    pointing staying within ``max_pointing_scatter`` degrees (which drops
    survey mosaics that reuse one object name).
    """
    names = ", ".join(_quote(i) for i in instruments)
    frames = query_frames(
        f"dp_cat='SCIENCE' AND dp_tech LIKE 'IMAGE%' AND instrument IN ({names}) "
        f"AND release_date < {_quote(released_before)}",
        columns=("instrument", "prog_id", "object", "filter_path", "mjd_obs", "ra", "dec", "exposure"))

    rows = []
    groups = frames.group_by(["instrument", "prog_id", "object", "filter_path", "night"])
    for key, group in zip(groups.groups.keys, groups.groups):
        t0, t1 = group["mjd_obs"].min(), group["mjd_obs"].max()
        if len(group) < min_frames or t1 - t0 < min_span_days:
            continue
        dec = float(np.median(group["dec"]))
        if (np.ptp(group["dec"]) > max_pointing_scatter
                or np.ptp(group["ra"]) * np.cos(np.radians(dec)) > max_pointing_scatter):
            continue
        rows.append((key["instrument"], key["prog_id"], key["object"], key["filter_path"],
                     int(key["night"]), len(group), t0, t1,
                     float(np.median(group["ra"])), dec, float(np.median(group["exposure"]))))
    return Table(rows=rows, names=("instrument", "prog_id", "object", "filter", "night", "n_frames",
                                   "mjd_start", "mjd_end", "ra", "dec", "exposure"))


def sequence_frames(instrument, object_name, night, filter_path=None, prog_id=None):
    """The SCIENCE frames of one sequence, in time order."""
    where = (f"instrument={_quote(instrument)} AND object={_quote(object_name)} AND dp_cat='SCIENCE' "
             f"AND mjd_obs >= {night + 0.5} AND mjd_obs < {night + 1.5}")
    if filter_path:
        where += f" AND filter_path={_quote(filter_path)}"
    if prog_id:
        where += f" AND prog_id={_quote(prog_id)}"
    table = query_frames(where)
    table.sort("mjd_obs")
    return table


def night_calibrations(instrument, night, wanted, window_days=1.0):
    """Raw calibration frames taken around a night.

    ``wanted`` maps an ESO ``dp_type`` (e.g. ``"BIAS"``, ``"FLAT,SKY"``) to a
    filter name, or to None for any filter (biases and darks don't depend on
    it). Frames within ``window_days`` either side of the night's midpoint are
    returned, closest in time first.
    """
    mid = night + 1.0
    table = query_frames(
        f"instrument={_quote(instrument)} AND dp_cat='CALIB' "
        f"AND mjd_obs >= {mid - window_days} AND mjd_obs < {mid + window_days}")
    keep = np.zeros(len(table), dtype=bool)
    for dp_type, filter_path in wanted.items():
        match = table["dp_type"] == dp_type
        if filter_path is not None:
            match &= table["filter_path"] == filter_path
        keep |= match
    table = table[keep]
    table = table[np.argsort(np.abs(table["mjd_obs"] - mid))]
    return table


def download(dp_ids, destination, workers=4, unzip=True):
    """Download frames into ``destination``, skipping any already there.

    Downloads run in ``workers`` threads, each with its own ``Eso`` client.
    Returns the local paths, in the order of ``dp_ids``.
    """
    from astroquery.eso import Eso

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)

    def local(dp_id):
        for suffix in (".fits", ".fits.Z", ".fits.gz", ".fits.fz"):
            path = destination / f"{dp_id}{suffix}"
            if path.exists():
                return path
        return None

    todo = [d for d in dp_ids if local(d) is None]
    log.info("%d of %d frames already present, %d to fetch", len(dp_ids) - len(todo), len(dp_ids), len(todo))

    def fetch(dp_id):
        try:
            Eso().retrieve_data([dp_id], destination=str(destination), unzip=unzip)
        except Exception:
            log.exception("failed to download %s", dp_id)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, _ in enumerate(pool.map(fetch, todo), 1):
            if i % 10 == 0 or i == len(todo):
                log.info("downloaded %d/%d", i, len(todo))

    paths = [local(d) for d in dp_ids]
    missing = [d for d, p in zip(dp_ids, paths) if p is None]
    if missing:
        log.warning("%d frames missing after download: %s", len(missing), missing[:5])
    return paths


def table_records(table, columns):
    """A table as a list of JSON-safe dicts."""
    out = []
    for row in table:
        record = {}
        for name in columns:
            value = row[name]
            if np.ma.is_masked(value):
                value = None
            elif isinstance(value, np.generic):
                value = value.item()
            record[name] = value
        out.append(record)
    return out


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(tmp, path)
