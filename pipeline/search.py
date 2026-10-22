"""Blind single-transit search over every light curve in a results JSON.

Usage
-----
    python -m pipeline.search RESULTS.json [--out CANDIDATES.json] [--plot-dir DIR] [--top 12]

Nothing here knows which star is the intended target: every source is
treated the same, and the target's rank is reported afterwards as a check.

For each star:

1. Cotrend: remove the strongest trends shared across the field (the
   leading principal components of the reference stars' normalised light
   curves: transparency, extinction, colour and position effects). A
   transit on one star isn't shared, so it survives.
2. Clip single-frame outliers (cosmic rays, satellites) against a running
   median.
3. Flag variable stars: smooth variability (robust standard deviation)
   far above the frame-to-frame scatter. They are listed separately, not
   ranked as transit candidates.
4. Slide a box of each trial duration across the night. The depth is the
   median out of the box minus the median inside it. It is scored against
   the noise on that timescale: the scatter of box-length averages taken
   outside the box, which includes correlated ("red") noise that
   frame-to-frame scatter misses. A single night can't constrain a period,
   so this is a single-event search.

Transit candidates are ranked by their best score. The scores are not calibrated false
alarm probabilities: the distribution of the top scores over the field is
what tells you how often noise alone reaches a given score.
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from . import store

log = logging.getLogger("pipeline.search")

DURATIONS_H = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
# Robust spread over frame-to-frame scatter above which a star counts as variable.
# White noise gives ~1; a 2.5% transit in 1% noise gives ~1.3.
VARIABLE_RATIO = 3.0
# Shared trends removed by cotrending.
N_COMPONENTS = 4
# Events within this time of each other (days) among the top candidates count as shared.
SHARED_EVENT_DAYS = 15.0 / 1440


def _array(values):
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def robust_p2p(x):
    """Point-to-point scatter, insensitive to slow trends and to outliers."""
    d = np.diff(x[np.isfinite(x)])
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)


def shared_trends(flux, usable, n_components=N_COMPONENTS):
    """Leading principal components (frames x n) of the ``usable`` stars' normalised light curves."""
    norm = flux[:, usable] / np.nanmedian(flux[:, usable], axis=0) - 1.0
    norm = norm / np.std(norm, axis=0)  # every reference star weighted equally
    u, sv, _ = np.linalg.svd(norm, full_matrices=False)
    explained = sv ** 2 / np.sum(sv ** 2)
    return u[:, :n_components], explained[:n_components]


def cotrend(f, basis):
    """``f`` (normalised to 1) with its least-squares fit to the shared trends removed."""
    design = np.column_stack([np.ones(f.size), basis])
    coeffs, *_ = np.linalg.lstsq(design, f, rcond=None)
    return f - basis @ coeffs[1:]


def clip_outliers(f, window=11, n_sigma=5.0):
    """``f`` with points more than ``n_sigma`` from a running median set to NaN."""
    from scipy.ndimage import median_filter

    smooth = median_filter(f, size=window, mode="nearest")
    resid = f - smooth
    sigma = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
    out = f.copy()
    out[np.abs(resid) > n_sigma * sigma] = np.nan
    return out


def timescale_noise(f, n):
    """Noise on an ``n``-point average: white plus correlated ("red") noise.

    The scatter of the means of every ``n``-point sliding window, never
    below the white-noise level sigma / sqrt(n). Sliding windows give a
    stable estimate even when only a few independent blocks fit.
    """
    f = f[np.isfinite(f)]
    if f.size < 2 * n:
        return np.nan
    white = robust_p2p(f) / np.sqrt(n)
    c = np.cumsum(np.insert(f, 0, 0.0))
    means = (c[n:] - c[:-n]) / n
    return float(max(np.std(means, ddof=1), white))


def box_search(t, f, durations_h=DURATIONS_H, step_h=0.05, min_points=5):
    """Best single box: (score, t_mid, duration_h, depth, n_in), scored against same-timescale noise."""
    best = (0.0, np.nan, np.nan, np.nan, 0)
    good = np.isfinite(f)
    t, f = t[good], f[good]
    for dur in durations_h:
        half = dur / 48.0
        for centre in np.arange(t.min() + half, t.max() - half + 1e-9, step_h / 24.0):
            inside = np.abs(t - centre) < half
            n_in = int(inside.sum())
            if n_in < min_points or n_in > 0.6 * t.size:
                continue
            depth = np.median(f[~inside]) - np.median(f[inside])
            noise = timescale_noise(f[~inside], n_in)
            if not np.isfinite(noise) or noise <= 0:
                continue
            score = depth / noise
            if score > best[0]:
                best = (float(score), float(centre), float(dur), float(depth), n_in)
    return best


def prepare_night(path):
    """Cotrended, clipped light curves of every star measured in every used frame of one night.

    Returns ``(t, rows)``: the used frames' BJD_TDB and, per star, a dict
    with its identity, ``detrended`` flux (normalised to 1, NaN where
    clipped), frame-to-frame ``scatter``, ``variability`` (robust spread
    over that scatter) and ``level`` (median brightness relative to the
    reference stars, for comparing nights).
    """
    data = store.read(path)
    frames, sources = data["frames"], data["sources"]
    t = _array(frames.get("bjd_tdb_date_obs", frames["bjd_tdb"]))
    used = np.array(frames["used"], dtype=bool)
    flux = np.stack([_array(s["flux"]) for s in sources], axis=1)[used]
    t = t[used]

    complete = np.all(np.isfinite(flux) & (flux > 0), axis=0)
    raw_scatter = np.array([robust_p2p(flux[:, j] / np.nanmedian(flux[:, j])) if complete[j] else np.nan
                            for j in range(flux.shape[1])])
    # Shared trends from the better half of the complete light curves.
    reference = complete & (raw_scatter <= np.nanmedian(raw_scatter[complete]))
    basis, explained = shared_trends(flux, reference)
    log.info("%d of %d stars measured in every frame; %d shared trends from %d reference stars "
             "explain %s of their variance", complete.sum(), len(sources), basis.shape[1], reference.sum(),
             ", ".join(f"{100 * e:.0f}%" for e in explained))

    # Brightness relative to the reference stars, for comparing nights (slow variables).
    ensemble = np.nansum(flux[:, reference], axis=1)
    rows = []
    for j, source in enumerate(sources):
        if not complete[j]:
            continue
        level = float(np.nanmedian(flux[:, j] / ensemble))
        f = clip_outliers(cotrend(flux[:, j] / np.nanmedian(flux[:, j]), basis))
        sigma = robust_p2p(f)
        spread = 1.4826 * np.nanmedian(np.abs(f - np.nanmedian(f)))
        rows.append({"index": j, "gaia_dr3": source["gaia_dr3"], "g_mag": source["g_mag"],
                     "role": source["role"], "ra": source["ra"], "dec": source["dec"],
                     "scatter": float(sigma), "variability": float(spread / sigma) if sigma > 0 else np.nan,
                     "level": level, "detrended": f})
    return t, rows


def search(path, top=12, plot_dir=None, out_path=None):
    t, rows = prepare_night(path)
    for r in rows:
        score, t_mid, dur, depth, n_in = box_search(t, r["detrended"])
        r.update(score=score, t_mid=t_mid, duration_h=dur, depth=depth, n_in=n_in)

    # A transit of a quiet star raises its spread only a little above its frame-to-frame
    # scatter; pulsators, contact binaries and slow variables raise it a lot.
    variables = sorted((r for r in rows if r["variability"] > VARIABLE_RATIO),
                       key=lambda r: r["variability"], reverse=True)
    rows = sorted((r for r in rows if not r["variability"] > VARIABLE_RATIO),
                  key=lambda r: r["score"], reverse=True)
    scores = np.array([r["score"] for r in rows])
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank
    for rank, r in enumerate(variables, 1):
        r["rank"] = f"V{rank}"
    # Events many stars share at the same time are systematics, not transits.
    times = np.array([r["t_mid"] for r in rows[:max(20, len(rows) // 20)]])
    for r in rows:
        r["shared_event"] = int(np.sum(np.abs(times - r["t_mid"]) < SHARED_EVENT_DAYS)) - 1 \
            if np.isfinite(r["t_mid"]) else 0
    target = next((r for r in rows + variables if r["role"] == "target"), None)
    log.info("%d variable stars (spread > %.1f x frame-to-frame scatter), %d searched for transits",
             len(variables), VARIABLE_RATIO, len(rows))
    for r in variables[:8]:
        log.info("variable %-4s G=%5.2f spread/scatter %5.1f  scatter %.2f%%  Gaia %d", r["rank"], r["g_mag"],
                 r["variability"], 100 * r["scatter"], r["gaia_dr3"])
    summary = {"n_searched": len(rows), "n_variable": len(variables),
               "variable_ratio_threshold": VARIABLE_RATIO,
               "score_percentiles": {p: float(np.percentile(scores, p)) for p in (50, 90, 99)}}
    if target:
        summary["target"] = {k: target.get(k) for k in ("rank", "score", "t_mid", "duration_h", "depth", "g_mag",
                                                         "scatter", "variability", "shared_event")}
        log.info("target ranks %d of %d (score %.1f; field median %.1f, 99th percentile %.1f)",
                 target["rank"], len(rows), target["score"], summary["score_percentiles"][50],
                 summary["score_percentiles"][99])
    for r in rows[:top]:
        log.info("#%-3d G=%5.2f %-10s score %6.1f  depth %6.2f%%  dur %.1f h  t_mid %.4f  scatter %.2f%%  "
                 "shared %2d  Gaia %d", r["rank"], r["g_mag"], r["role"], r["score"], 100 * r["depth"],
                 r["duration_h"], r["t_mid"], 100 * r["scatter"], r["shared_event"], r["gaia_dr3"])

    if out_path:
        strip = lambda rs: [{k: v for k, v in r.items() if k != "detrended"} for r in rs]
        Path(out_path).write_text(json.dumps({"results": str(path), "summary": summary,
                                              "candidates": strip(rows), "variables": strip(variables)},
                                             indent=1) + "\n")
    if plot_dir:
        _plot(rows[:top], t, Path(plot_dir), target)
        _plot(variables[:top], t, Path(plot_dir), None, name="variables.png")
    return rows, summary


def _plot(rows, t, plot_dir, target, name="candidates.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir.mkdir(parents=True, exist_ok=True)
    shown = list(rows)
    if target and target not in shown:
        shown.append(target)
    n = len(shown)
    cols = 3
    fig, axes = plt.subplots((n + cols - 1) // cols, cols, figsize=(12, 2.4 * ((n + cols - 1) // cols)),
                             sharex=True, squeeze=False)
    hours = (t - t.min()) * 24
    for ax, r in zip(axes.ravel(), shown):
        ax.plot(hours, r["detrended"], ".", ms=2, color="0.4")
        if np.isfinite(r["t_mid"]):
            c = (r["t_mid"] - t.min()) * 24
            ax.axvspan(c - r["duration_h"] / 2, c + r["duration_h"] / 2, color="C3", alpha=0.15)
        label = " (target)" if r["role"] == "target" else ""
        ax.set_title(f"#{r['rank']} G={r['g_mag']:.1f}{label}: score {r['score']:.1f}, "
                     f"depth {100 * r['depth']:.1f}%, var {r['variability']:.1f}", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.supxlabel("hours from first frame", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_dir / name, dpi=110)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results")
    parser.add_argument("--out")
    parser.add_argument("--plot-dir")
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    search(args.results, top=args.top, plot_dir=args.plot_dir, out_path=args.out)


if __name__ == "__main__":
    main()
