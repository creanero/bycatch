"""JSON results store: one file per extracted sequence (per chip, for mosaics).

Layout of a results file::

    {
      "schema": "bycatch-lightcurves/1",
      "sequence": {...},     # archive, instrument, target, night, site, files, apertures,
                             # gain/read noise, astrometry, software versions
      "frames": {"bjd_tdb": [...], "jd_utc": [...], "airmass": [...], "aligned": [...]},
      "sources": [           # one per photometered star, target first
        {"gaia_dr3": 123, "g_mag": 12.3, "ra": ..., "dec": ..., "role": "target" | "comparison" | "field",
         "flux": [...], "uncertainty": [...], "x": [...], "y": [...]},
        ...
      ],
      "analysis": {...}      # e.g. comparison-star choice and the transit fit
    }

NaN is written as null. The fields map onto the kdb+ schema in database/:
``sources`` (gaia_dr3, ra, dec), ``obs_sessions`` (sequence) and
``observations`` (frames x sources).
"""

import json
import math
import os
from pathlib import Path

import numpy as np

SCHEMA = "bycatch-lightcurves/1"


def _clean(value):
    if isinstance(value, dict):
        return {k: _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return _clean(value.tolist())
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write(path, sequence, frames, sources, analysis=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {"schema": SCHEMA, "sequence": sequence, "frames": frames, "sources": sources,
            "analysis": analysis or {}}
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(_clean(data), allow_nan=False) + "\n")
    os.replace(tmp, path)
    return path


def read(path):
    data = json.loads(Path(path).read_text())
    if data.get("schema") != SCHEMA:
        raise ValueError(f"{path}: unknown schema {data.get('schema')!r}")
    return data
