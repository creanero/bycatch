"""Resolve a target definition to archive frames and download them.

Usage
-----
    python -m access.fetch targets/pots1_wfi.json [--data-root ~/bycatch-data]
        [--manifest-only] [--workers 4]

Writes ``<data-root>/<archive>/<instrument>/<night>_<object>/manifest.json``
listing every science and calibration frame (with its ESO metadata), then
downloads the frames into ``science/`` and ``calib/`` beside it. Re-running
skips frames already downloaded, so an interrupted fetch can just be
restarted.
"""

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from . import eso

MANIFEST_COLUMNS = ("dp_id", "dp_cat", "dp_type", "filter_path", "exposure", "mjd_obs", "ra", "dec",
                    "access_estsize", "release_date")


def sequence_dir(data_root, target):
    safe_object = target["object"].replace("/", "_").replace(" ", "_")
    return (Path(data_root).expanduser() / target["archive"].lower() / target["instrument"].lower()
            / f"{target['night']}_{safe_object}")


def resolve(target):
    """The manifest for a target definition: its science and calibration frames."""
    if target["archive"] != "ESO":
        raise ValueError(f"unsupported archive {target['archive']!r}")
    science = eso.sequence_frames(target["instrument"], target["object"], target["night"],
                                  filter_path=target.get("filter"), prog_id=target.get("prog_id"))
    wanted = {dp_type: (target.get("filter") if same_filter else None)
              for dp_type, same_filter in target.get("calibrations", {}).items()}
    calib = eso.night_calibrations(target["instrument"], target["night"], wanted)
    return {
        "target": target,
        "resolved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "science": eso.table_records(science, MANIFEST_COLUMNS),
        "calibrations": eso.table_records(calib, MANIFEST_COLUMNS),
    }


def summarise(manifest):
    science, calib = manifest["science"], manifest["calibrations"]
    size_gb = sum((f["access_estsize"] or 0) for f in science + calib) / 1024 ** 2
    counts = {}
    for f in calib:
        counts[f"{f['dp_type']} {f['filter_path']}"] = counts.get(f"{f['dp_type']} {f['filter_path']}", 0) + 1
    span_h = (science[-1]["mjd_obs"] - science[0]["mjd_obs"]) * 24 if science else 0.0
    logging.info("%d science frames over %.2f h; calibrations %s; ~%.1f GB compressed",
                 len(science), span_h, counts, size_gb)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("target", help="target definition JSON")
    parser.add_argument("--data-root", default="~/bycatch-data")
    parser.add_argument("--manifest-only", action="store_true", help="resolve and write the manifest, no download")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    target = json.loads(Path(args.target).read_text())
    out = sequence_dir(args.data_root, target)
    manifest = resolve(target)
    eso.write_json(out / "manifest.json", manifest)
    logging.info("manifest: %s", out / "manifest.json")
    summarise(manifest)
    if args.manifest_only:
        return

    eso.download([f["dp_id"] for f in manifest["calibrations"]], out / "calib", workers=args.workers)
    eso.download([f["dp_id"] for f in manifest["science"]], out / "science", workers=args.workers)
    logging.info("done: %s", out)


if __name__ == "__main__":
    main()
