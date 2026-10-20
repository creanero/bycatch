"""Refit the target's transit in a results JSON, without re-extracting.

Usage
-----
    python -m pipeline.refit RESULTS.json [--plot OUT.png]

Uses the stored relative light curve, frame times, airmass and planet
parameters, refits with ``lightcurve.fit_with_airmass``, and rewrites
``analysis.transit_fit`` in place. MicroObservatory results use the
DATE-OBS-based times when present.
"""

import argparse
import logging

import numpy as np

from . import store
from .lightcurve import LIMB_DARKENING, fit_with_airmass

log = logging.getLogger("pipeline.refit")


def _array(values):
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def refit(path, plot_path=None):
    data = store.read(path)
    frames, analysis = data["frames"], data["analysis"]
    planet = analysis["planet"]
    times = _array(frames.get("bjd_tdb_date_obs", frames["bjd_tdb"]))
    used = np.array(frames["used"], dtype=bool)
    relative = _array(analysis["relative_flux"])
    airmass = _array(frames["airmass"])

    values, errors, extra = fit_with_airmass(times[used], relative[used], airmass[used], planet)
    fit = {"tmid_bjd": values["tmid"], "tmid_err_min": errors["tmid"] * 1440, "rprs": values["rprs"],
           "rprs_err": errors["rprs"], "scale": values["scale"],
           "o_minus_c_min": (values["tmid"] - planet["tmid_predicted"]) * 1440,
           "airmass_poly": extra["airmass_poly"], "residual_rms": extra["rms"],
           "n_points": int(extra["kept"].sum()), "n_clipped": extra["n_clipped"],
           "method": "fit_with_airmass: relative times, grid start, 4-sigma clipping"}
    analysis["transit_fit"] = fit
    store.write(path, data["sequence"], frames, data["sources"], analysis)
    log.info("%s: T_mid %.5f +- %.2f min, O-C %+.2f min, Rp/Rs %.3f +- %.3f, rms %.4f, %d clipped",
             planet["name"], fit["tmid_bjd"], fit["tmid_err_min"], fit["o_minus_c_min"], fit["rprs"],
             fit["rprs_err"], fit["residual_rms"], fit["n_clipped"])

    if plot_path:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = (times[used] - planet["tmid_predicted"]) * 24
        order = np.argsort(t)
        kept = extra["kept"]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(t[kept], extra["detrended"][kept], ".", color="0.45", ms=4, label="relative flux, airmass-detrended")
        if (~kept).any():
            ax.plot(t[~kept], extra["detrended"][~kept], "x", color="0.7", ms=4, label="clipped (>4σ)")
        tk, fk = t[kept], extra["detrended"][kept]
        o = np.argsort(tk)
        nbin = 5 if len(tk) > 100 else 1
        if nbin > 1:
            nb = len(tk) // nbin * nbin
            ax.plot(tk[o][:nb].reshape(-1, nbin).mean(1), fk[o][:nb].reshape(-1, nbin).mean(1), "o",
                    color="C0", ms=5, label=f"binned ×{nbin}")
        ax.plot(t[order], extra["model"][order], "-", color="C3", label="transit fit")
        ax.axvline(0, color="0.6", ls="--", lw=1, label="predicted $T_{mid}$")
        ax.axvline(fit["o_minus_c_min"] / 60, color="C3", ls=":", lw=1, label="fitted $T_{mid}$")
        ax.set_xlabel("hours from predicted mid-transit")
        ax.set_ylabel("relative flux")
        ax.set_title(f"{planet['name']}: $T_{{mid}}$ O−C {fit['o_minus_c_min']:+.1f} ± {fit['tmid_err_min']:.1f} min, "
                     f"$R_p/R_s$ {fit['rprs']:.3f} ± {fit['rprs_err']:.3f}", fontsize=10)
        ax.legend(fontsize=7, loc="lower left")
        fig.tight_layout()
        fig.savefig(plot_path, dpi=130)
    return fit


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results")
    parser.add_argument("--plot")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    refit(args.results, args.plot)


if __name__ == "__main__":
    main()
