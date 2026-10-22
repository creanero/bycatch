"""Collate per-chip period-search results into one ranked list for the whole field.

Usage
-----
    python -m pipeline.collate SEARCH_DIR [SEARCH_DIR ...] --out FIELD.json [--top 25]

Each ``SEARCH_DIR`` holds the ``candidates.json`` written by
``pipeline.periodsearch`` for one chip. Candidates are ranked across all
chips the same way (strictly vetted first, then vetted, then by SNR), and
each candidate and excluded variable is matched against Gaia DR3's
variability classification (``gaiadr3.vari_classifier_result``), so stars
Gaia already knows as eclipsing binaries or other variables are labelled.
"""

import argparse
import json
import logging
import time
from pathlib import Path

log = logging.getLogger("pipeline.collate")


def gaia_variability(source_ids):
    """{source_id: (best_class_name, best_class_score)} from Gaia DR3 for the given IDs."""
    if not source_ids:
        return {}
    from astroquery.gaia import Gaia

    out = {}
    ids = sorted(set(int(s) for s in source_ids))
    for start in range(0, len(ids), 500):
        chunk = ids[start:start + 500]
        query = ("SELECT source_id, best_class_name, best_class_score FROM gaiadr3.vari_classifier_result "
                 f"WHERE source_id IN ({', '.join(map(str, chunk))})")
        for attempt in range(5):
            try:
                rows = Gaia.launch_job_async(query).get_results()
                break
            except Exception as exc:
                if attempt == 4:
                    raise
                log.warning("Gaia query failed (%s); retrying", exc)
                time.sleep(30 * (attempt + 1))
        for r in rows:
            out[int(r["source_id"])] = (str(r["best_class_name"]), float(r["best_class_score"]))
    return out


def collate(search_dirs, out_path, top=25):
    candidates, variables, summaries = [], [], {}
    for d in search_dirs:
        path = Path(d) / "candidates.json"
        if not path.exists():
            log.warning("no results in %s", d)
            continue
        data = json.loads(path.read_text())
        chip = Path(d).name
        summaries[chip] = data["summary"]
        for c in data["candidates"]:
            candidates.append(dict(c, chip=chip))
        variables += [{"gaia_dr3": g, "chip": chip} for g in data.get("variable_excluded", [])]

    known = gaia_variability([c["gaia_dr3"] for c in candidates] + [v["gaia_dr3"] for v in variables])
    for c in candidates + variables:
        c["gaia_variability"] = known.get(int(c["gaia_dr3"]))
    candidates.sort(key=lambda c: (not c.get("vetted_strict"), not c.get("vetted"), -c["snr"]))
    for rank, c in enumerate(candidates, 1):
        c["field_rank"] = rank

    n_strict = sum(bool(c.get("vetted_strict")) for c in candidates)
    log.info("%d chips, %d multi-night candidates (%d strictly vetted), %d variables excluded; "
             "%d candidates and %d variables have a Gaia DR3 variability class",
             len(summaries), len(candidates), n_strict, len(variables),
             sum(c["gaia_variability"] is not None for c in candidates),
             sum(v["gaia_variability"] is not None for v in variables))
    for c in candidates[:top]:
        vari = f"{c['gaia_variability'][0]} ({c['gaia_variability'][1]:.2f})" if c["gaia_variability"] else "-"
        log.info("#%-3d %-12s %s G=%5.2f %-6s SNR %5.1f  P %8.5f d  depth %5.2f%%  dur %.1f h  Gaia var %s  %d",
                 c["field_rank"], c["chip"], "S" if c.get("vetted_strict") else ("V" if c.get("vetted") else " "),
                 c["g_mag"], c["role"], c["snr"], c["period"], 100 * c["depth"], c["duration_h"], vari,
                 c["gaia_dr3"])
    for v in variables:
        if v["gaia_variability"]:
            log.info("excluded variable %s %d: Gaia %s (%.2f)", v["chip"], v["gaia_dr3"], *v["gaia_variability"])

    Path(out_path).write_text(json.dumps({"chips": summaries, "candidates": candidates, "variables": variables},
                                         indent=1, default=float) + "\n")
    return candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("search_dirs", nargs="+")
    parser.add_argument("--out", required=True)
    parser.add_argument("--top", type=int, default=25)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    collate(args.search_dirs, args.out, top=args.top)


if __name__ == "__main__":
    main()
