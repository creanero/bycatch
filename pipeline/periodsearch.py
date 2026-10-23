"""Blind periodic transit search across several nights of the same field.

Usage
-----
    python -m pipeline.periodsearch NIGHT1.json NIGHT2.json [...] [--out CANDIDATES.json]
        [--plot-dir DIR] [--min-period 0.5] [--max-period 10] [--workers 16]

Stars are matched across nights by Gaia DR3 ID. Each night is prepared as
in ``pipeline.search`` (cotrended, clipped, variables flagged), and stars
that are variable on any night are left out. Each remaining star's nights
are joined and searched with a box least squares periodogram (astropy).
A real transiting planet dips again at the same depth and duration one or
more orbits later; a one-night systematic doesn't repeat in step.

Per-point uncertainties are each star's frame-to-frame scatter on each
night, inflated by that night's correlated-noise factor at a 1.5 h
timescale. The factor is the median over the field's stars, not the star's
own: a star's own factor is inflated by its transit, if it has one.

Trial periods within 2% of a whole number of days, or of 1/2 or 1/3 day,
are masked. A systematic at the same time of night on every night
(extinction, flexure, scattered light) lines up at those periods.

Nothing here knows which star is the intended target. Stars are ranked by
the periodogram's best signal-to-noise, and only candidates with in-transit
points on at least two nights count. Each candidate is then vetted, as you
would by eye:

- at least two events are significant on their own (``EVENT_MIN_SNR``),
  and no predicted event shows the star clearly not dimming;
- the events' depths agree within ``DEPTH_AGREEMENT_SIGMA``;
- at least one event is fully inside its night, with at least
  ``BASELINE_H`` of data before and after it, not cut off at the start or
  end where trends dominate (``vetted``); at least two such complete
  events make it ``vetted_strict``.

With many nights a real planet has some partial or noisy events, so the
checks count events rather than requiring every one to pass. With two
nights the two readings are the same.

Stars whose brightness relative to the field changes between nights by
more than ``LEVEL_CHANGE_SIGMA`` times the field's spread (slow variables,
such as long-period variables, which look flat within a night) are
excluded with the variables.

Candidates passing all three are listed first (``vetted``). With a few nights the period is
ambiguous: the events fix the spacing between transits, and any whole
fraction of it in the search range fits. The aliases are listed with each
candidate.
"""

import argparse
import json
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.timeseries import BoxLeastSquares

from .search import box_search, prepare_night, robust_p2p, timescale_noise

log = logging.getLogger("pipeline.periodsearch")

DURATIONS_H = (1.0, 1.5, 2.0, 3.0)
DAILY_ALIAS_TOLERANCE = 0.02
# Smooth variability left after removing each night's best box, on any night.
MULTI_NIGHT_VARIABLE_RATIO = 2.0
EVENT_MIN_SNR = 3.0
BASELINE_H = 0.5
LEVEL_CHANGE_SIGMA = 5.0
LEVEL_CHANGE_MIN = 0.02
DEPTH_AGREEMENT_SIGMA = 3.0
RED_NOISE_TIMESCALE_H = 1.5


def variability_beyond_box(t, star):
    """Robust spread over frame-to-frame scatter once the night's best box is taken out.

    A transit raises a light curve's spread, sometimes a lot when it covers
    much of a short night, but one box removes it. A variable star's spread
    survives.
    """
    f = star["detrended"].copy()
    _, t_mid, dur, depth, _ = box_search(t, f)
    if np.isfinite(t_mid):
        f[np.abs(t - t_mid) < dur / 48] += depth
    sigma = robust_p2p(f)
    spread = 1.4826 * np.nanmedian(np.abs(f - np.nanmedian(f)))
    return spread / sigma if sigma > 0 else np.inf


def level_changes(stars, by_gaia):
    """Stars whose brightness relative to the field changes between nights more than the field allows."""
    stars = sorted(stars)
    if len(by_gaia) < 2:
        return set()
    levels = np.array([[night[g]["level"] for g in stars] for night in by_gaia])  # nights x stars
    ratio = np.log(levels / levels[0])
    ratio -= np.median(ratio, axis=1, keepdims=True)  # nights differ in aperture, seeing, transparency
    spread = 1.4826 * np.median(np.abs(ratio[1:]), axis=1, keepdims=True)
    limit = np.maximum(LEVEL_CHANGE_SIGMA * spread, LEVEL_CHANGE_MIN)
    changed = np.any(np.abs(ratio[1:]) > limit, axis=0)
    return {g for g, c in zip(stars, changed) if c}


def red_noise_factor(t, f, timescale_h=RED_NOISE_TIMESCALE_H):
    """How much worse than white the noise is on ``timescale_h``: (block-average scatter) / (sigma / sqrt(n))."""
    good = np.isfinite(f)
    cadence = np.median(np.diff(t[good]))
    n = max(2, int(round(timescale_h / 24 / cadence)))
    sigma = robust_p2p(f)
    noise = timescale_noise(f, n)
    return max(1.0, noise / (sigma / np.sqrt(n))) if np.isfinite(noise) and sigma > 0 else 1.0


_JOB = {}


def _search_star(item):
    gaia, t, f, dy = item
    model = BoxLeastSquares(t, f, dy)
    durations = np.array(DURATIONS_H) / 24
    result = model.power(_JOB["periods"], durations, objective="snr")
    i = int(np.argmax(result.power))
    period, t0, duration = result.period[i], result.transit_time[i], result.duration[i]
    stats = model.compute_stats(period, duration, t0)
    nights_hit = len({int(np.floor(tt - 0.5)) for tt, n in zip(stats["transit_times"], stats["per_transit_count"])
                      if n > 0})
    return {"gaia_dr3": gaia, "snr": float(result.power[i]), "period": float(period), "t0": float(t0),
            "duration_h": float(duration * 24), "depth": float(result.depth[i]),
            "depth_err": float(result.depth_err[i]), "nights_with_transit": nights_hit,
            "depth_odd": float(stats["depth_odd"][0]), "depth_even": float(stats["depth_even"][0]),
            "transit_times": [float(x) for x, n in zip(stats["transit_times"], stats["per_transit_count"]) if n > 0]}


def daily_alias(periods, tolerance=DAILY_ALIAS_TOLERANCE):
    """True where a period is within ``tolerance`` (fractional) of n days or 1/n day (n = 1, 2, 3, ...)."""
    near_days = np.abs(periods / np.maximum(np.round(periods), 1) - 1) < tolerance
    near_fraction = (periods < 1) & (np.abs(1 / periods / np.round(1 / periods) - 1) < tolerance)
    return (near_days & (periods >= 0.98)) | near_fraction


def vet(r, nights, by_gaia, betas):
    """Per-night events of a star's best solution, and whether they pass the vetting checks."""
    events = []
    half = r["duration_h"] / 48
    for (t, _), night, beta in zip(nights, by_gaia, betas):
        f = night[r["gaia_dr3"]]["detrended"]
        good = np.isfinite(f)
        t, f = t[good], f[good]
        n = np.arange(np.ceil((t.min() - r["t0"]) / r["period"] - 1), np.floor((t.max() - r["t0"]) / r["period"]) + 2)
        for tc in r["t0"] + n * r["period"]:
            inside = np.abs(t - tc) < half
            if inside.sum() < 3:
                continue
            depth = float(np.median(f[~inside]) - np.mean(f[inside]))
            err = float(night[r["gaia_dr3"]]["scatter"] * beta / np.sqrt(inside.sum()))
            edge = bool(tc - half - BASELINE_H / 24 < t.min() or tc + half + BASELINE_H / 24 > t.max())
            events.append({"t_mid": float(tc), "n_in": int(inside.sum()), "depth": depth, "depth_err": err,
                           "edge": edge})
    if len(events) < 2:
        return {"events": events, "vetted": False, "vetted_strict": False}
    depths = np.array([e["depth"] for e in events])
    errs = np.array([e["depth_err"] for e in events])
    mean = np.sum(depths / errs ** 2) / np.sum(1 / errs ** 2)
    consistent = bool(np.all(np.abs(depths - mean) < DEPTH_AGREEMENT_SIGMA * errs))
    snr = depths / errs
    significant = bool(np.sum(snr > EVENT_MIN_SNR) >= 2 and not np.any(snr < -EVENT_MIN_SNR))
    n_complete = sum(not e["edge"] for e in events)
    vetted = significant and consistent and n_complete >= 1
    return {"events": events, "events_significant": significant, "depths_consistent": consistent,
            "n_complete_events": n_complete, "vetted": vetted, "vetted_strict": vetted and n_complete >= 2}


def aliases(transit_times, min_period, max_period):
    """Periods in range consistent with the observed transit times (gap / N for whole N)."""
    if len(transit_times) < 2:
        return []
    gap = transit_times[-1] - transit_times[0]
    n = np.arange(max(1, int(np.ceil(gap / max_period))), int(np.floor(gap / min_period)) + 1)
    return [float(gap / k) for k in n]


def period_search(paths, out_path=None, plot_dir=None, min_period=0.7, max_period=10.0, workers=16, top=12):
    nights = [prepare_night(p) for p in paths]
    by_gaia = [{r["gaia_dr3"]: r for r in rows} for _, rows in nights]
    common = set(by_gaia[0]).intersection(*by_gaia[1:])
    variable = {g for g in common
                if any(variability_beyond_box(t, night[g]) > MULTI_NIGHT_VARIABLE_RATIO
                       for (t, _), night in zip(nights, by_gaia))}
    slow = level_changes(common, by_gaia)
    log.info("%d stars change brightness between nights (slow variables)", len(slow - variable))
    variable |= slow
    stars = sorted(common - variable)
    log.info("%d stars in all %d nights; %d variable on some night, %d searched", len(common), len(paths),
             len(variable), len(stars))

    betas = [float(np.median([red_noise_factor(t, night[g]["detrended"]) for g in stars]))
             for (t, _), night in zip(nights, by_gaia)]
    log.info("field correlated-noise factor at %.1f h per night: %s", RED_NOISE_TIMESCALE_H,
             ", ".join(f"{b:.2f}" for b in betas))
    items = []
    for gaia in stars:
        ts, fs, dys = [], [], []
        for (t, _), night, beta in zip(nights, by_gaia, betas):
            f = night[gaia]["detrended"]
            good = np.isfinite(f)
            ts.append(t[good])
            fs.append(f[good])
            dys.append(np.full(good.sum(), night[gaia]["scatter"] * beta))
        items.append((gaia, np.concatenate(ts), np.concatenate(fs), np.concatenate(dys)))

    # Every star shares the frame times (to within clipping), so one period grid serves all.
    durations = np.array(DURATIONS_H) / 24
    periods = BoxLeastSquares(items[0][1], items[0][2]).autoperiod(
        durations, minimum_period=min_period, maximum_period=max_period, frequency_factor=1.0)
    masked = daily_alias(periods)
    periods = periods[~masked]
    _JOB["periods"] = periods
    log.info("%d trial periods %.2f-%.1f d (%d near daily aliases masked), durations %s h, %d workers",
             periods.size, periods.min(), periods.max(), masked.sum(), DURATIONS_H, workers)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_search_star, items, chunksize=8))

    info = by_gaia[0]
    for r in results:
        r.update(vet(r, nights, by_gaia, betas))
        r.update({k: info[r["gaia_dr3"]][k] for k in ("g_mag", "role", "ra", "dec")})
        r["scatter"] = [night[r["gaia_dr3"]]["scatter"] for night in by_gaia]
        r["aliases"] = aliases(r["transit_times"], min_period, max_period)
    multi = sorted((r for r in results if r["nights_with_transit"] >= 2),
                   key=lambda r: (not r["vetted_strict"], not r["vetted"], -r["snr"]))
    single = [r for r in results if r["nights_with_transit"] < 2]
    for rank, r in enumerate(multi, 1):
        r["rank"] = rank
    snrs = np.array([r["snr"] for r in multi])
    summary = {"nights": [str(p) for p in paths], "n_searched": len(results), "n_multi_night": len(multi),
               "n_variable_excluded": len(variable), "red_noise_factor": betas,
               "period_range": [min_period, max_period], "daily_alias_tolerance": DAILY_ALIAS_TOLERANCE,
               "snr_percentiles": {p: float(np.percentile(snrs, p)) for p in (50, 90, 99)} if len(snrs) else {}}
    summary["n_vetted"] = sum(r["vetted"] for r in multi)
    summary["n_vetted_strict"] = sum(r["vetted_strict"] for r in multi)
    log.info("%d of %d best solutions put transits on 2+ nights (%d pass vetting, %d with 2+ complete events); "
             "SNR median %.1f, 99th percentile %.1f", len(multi), len(results), summary["n_vetted"],
             summary["n_vetted_strict"],
             summary["snr_percentiles"].get(50, np.nan), summary["snr_percentiles"].get(99, np.nan))
    target = next((r for r in results if r["role"] == "target"), None)
    if target:
        summary["target"] = {k: target.get(k) for k in ("rank", "snr", "period", "t0", "duration_h", "depth",
                                                         "nights_with_transit", "transit_times", "vetted",
                                                         "vetted_strict", "events")}
        log.info("target: rank %s, SNR %.1f, P %.5f d, depth %.2f%%, transits on %d nights, vetted %s, strict %s",
                 target.get("rank", "-"), target["snr"], target["period"], 100 * target["depth"],
                 target["nights_with_transit"], target["vetted"], target["vetted_strict"])
    for r in multi[:top]:
        log.info("#%-3d %s G=%5.2f %-6s SNR %6.1f  P %8.5f d  depth %5.2f%%  dur %.1f h  events %s  Gaia %d",
                 r["rank"], "S" if r["vetted_strict"] else ("V" if r["vetted"] else " "), r["g_mag"], r["role"], r["snr"], r["period"],
                 100 * r["depth"], r["duration_h"],
                 " ".join(f"{100 * e['depth']:.2f}±{100 * e['depth_err']:.2f}%{'(edge)' if e['edge'] else ''}"
                          for e in r["events"]), r["gaia_dr3"])

    if out_path:
        Path(out_path).write_text(json.dumps({"summary": summary, "candidates": multi,
                                              "single_night": single, "variable_excluded": sorted(variable)},
                                             indent=1, default=float) + "\n")
    if plot_dir:
        _plot(multi[:top], target, nights, by_gaia, Path(plot_dir))
    return multi, summary


def _plot(rows, target, nights, by_gaia, plot_dir):
    """One row per candidate: each night's light curve, with the predicted transit windows shaded."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir.mkdir(parents=True, exist_ok=True)
    shown = list(rows) + ([target] if target and target not in rows else [])
    fig, axes = plt.subplots(len(shown), len(nights), figsize=(4.2 * len(nights), 1.7 * len(shown)),
                             squeeze=False, sharey="row")
    for i, r in enumerate(shown):
        for k, ((t, _), night) in enumerate(zip(nights, by_gaia)):
            ax = axes[i, k]
            f = night[r["gaia_dr3"]]["detrended"]
            ax.plot((t - t.min()) * 24, f, ".", ms=1.5, color="0.4")
            n = np.arange(np.ceil((t.min() - r["t0"]) / r["period"] - 1), np.floor((t.max() - r["t0"]) / r["period"]) + 2)
            for tc in r["t0"] + n * r["period"]:
                c = (tc - t.min()) * 24
                ax.axvspan(c - r["duration_h"] / 2, c + r["duration_h"] / 2, color="C3", alpha=0.15)
            ax.set_xlim(0, (t.max() - t.min()) * 24)
            ax.tick_params(labelsize=6)
            if k == 0:
                label = (" target" if r["role"] == "target" else "") + \
                    (" ✓✓" if r.get("vetted_strict") else " ✓" if r.get("vetted") else "")
                ax.set_ylabel(f"#{r.get('rank', '-')}{label}\nG={r['g_mag']:.1f}", fontsize=7)
            if i == 0:
                ax.set_title(f"night {int(np.floor(t.min() - 2400000.5 - 0.5))}", fontsize=8)
        axes[i, -1].text(1.02, 0.5, f"SNR {r['snr']:.1f}\nP {r['period']:.4f} d\n{100 * r['depth']:.1f}%",
                         transform=axes[i, -1].transAxes, fontsize=7, va="center")
    fig.supxlabel("hours from first frame of each night", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_dir / "period_candidates.png", dpi=110)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="+")
    parser.add_argument("--out")
    parser.add_argument("--plot-dir")
    parser.add_argument("--min-period", type=float, default=0.7)
    parser.add_argument("--max-period", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    period_search(args.results, out_path=args.out, plot_dir=args.plot_dir, min_period=args.min_period,
                  max_period=args.max_period, workers=args.workers, top=args.top)


if __name__ == "__main__":
    main()
