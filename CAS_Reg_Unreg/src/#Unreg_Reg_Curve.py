#Unreg_Reg_Curve.py
# -*- coding: utf-8 -*-
"""
Unregulated-regulated scatter and the inferred regulated frequency curve.

PART 1 -- SCATTER
    Unregulated peak against adjusted regulated peak at Castle Rock, drawn twice
    (arithmetic and log-log), both with a 1:1 line and the largest events called
    out by water year. The 1:1 line is the no-reservoir reference: points below
    it are the reduction the project achieved.

PART 2 -- INFERRED REGULATED FREQUENCY CURVE
    The regulated curve is NOT fitted analytically. Regulated peaks do not
    follow an analytical distribution -- operating rules put hard breaks in
    them -- so the AEP is inherited from the unregulated side instead:

        unregulated AEP  ->  unregulated peak (Expected curve)
                         ->  unreg-reg relationship
                         ->  regulated peak at that AEP

    Peak-to-peak is used because peak and 1-day tied as the critical duration
    (log r-squared 0.850 vs 0.836, not distinguishable at n=44), and a
    peak-to-peak transform is the easier one to explain.

    The transform line is the log-log power law fitted to the observed points.
    It is a ROUGH line through the scatter, plotted with a band at +/- one
    standard error, and is meant to be replaced by a hand-drawn curve once the
    synthetics populate the upper end.

    The frequency plot is formatted in the HEC-SSP idiom: normal-probability
    x-axis increasing to the right in AEP, log flow on y, both records plotted
    as points at their plotting positions.

CAUTION
    The transform is fitted over unregulated 22,000-218,000 cfs. Beyond that it
    is extrapolation, and the fitted exponent below 1 makes the curve flatten,
    so the extrapolated regulated flows look increasingly favourable. The plot
    marks the fitted range explicitly.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
DATASET_CSV = r"../output/diagnostics/critical_duration_adjusted_dataset.csv"
FITS_CSV = r"../output/critical_duration_adjusted_fits.csv"
UNREG_FREQ_CSV = r"../../CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv"

OUT_DIR = r"../output"
PLOT_STEM = r"../output/diagnostics/unreg_reg"

UNREG_COL = "unreg_Peak_1-hr"      # peak-to-peak transform
FREQ_DURATION = "Peak"
FREQ_VALUE_COL = "Expected"        # the expected-probability curve

CALLOUT_TOP_N = 8                  # label this many largest events
CALLOUT_MIN_CFS = 90000.0          # ...and anything above this

# How the observed points are placed on the frequency plot.
#   "from_curve" : each year's AEP is read off the UNREGULATED frequency curve
#                  at its own unregulated peak, and the regulated peak is drawn
#                  at that same AEP. This is the method itself made visible --
#                  the regulated point inherits the unregulated AEP -- and it is
#                  the only consistent choice here, because these 44 years are
#                  the regulated-era subset of a 98-year record. Weibull
#                  positions computed on 44 years would place them at far higher
#                  AEPs than the curve, which is fitted to all 98.
#   "weibull"    : plotting positions from this 44-year sample alone. Kept for
#                  comparison; expect the points to sit left of the curve.
PLOTTING_BASIS = "from_curve"

# SSP-style frequency axis
AEP_TICKS = [0.999, 0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02,
             0.01, 0.005, 0.002, 0.001]
AEP_LIMITS = (0.999, 0.0005)
FLOW_LIMITS = (10000.0, 400000.0)

# ----------------------------------------------------------------------------


def weibull_plotting_positions(values):
    """Weibull plotting positions (i/(n+1)) as AEP, largest first."""
    clean = np.sort(np.asarray(values, dtype=float))[::-1]
    n = len(clean)
    return clean, (np.arange(1, n + 1)) / float(n + 1)


def fit_power_law(x, y):
    """log10(y) on log10(x). Returns a, b, r2 and the log-space standard error."""
    good = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    lx, ly = np.log10(x[good]), np.log10(y[good])
    fit = stats.linregress(lx, ly)
    resid = ly - (fit.slope * lx + fit.intercept)
    return {"a": 10 ** fit.intercept, "b": fit.slope, "r2": fit.rvalue ** 2,
            "se_dex": float(np.std(resid, ddof=2)), "n": int(good.sum()),
            "x_min": float(x[good].min()), "x_max": float(x[good].max())}


def apply_power_law(fit, x):
    return fit["a"] * np.asarray(x, dtype=float) ** fit["b"]


def annotate_points(ax, x, y, labels, log_axes):
    """Callouts placed alternately above and below to reduce collisions."""
    for i, (xi, yi, text) in enumerate(zip(x, y, labels)):
        dy = 14 if i % 2 == 0 else -18
        ax.annotate(text, (xi, yi), xytext=(10, dy), textcoords="offset points",
                    fontsize=8, color="#c0392b",
                    arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.7))


def plot_scatter(data, fit, stem):
    """Arithmetic and log-log scatter, each with a 1:1 line and callouts."""
    x = data[UNREG_COL].values
    y = data["reg_peak"].values
    big = data.nlargest(CALLOUT_TOP_N, UNREG_COL)
    big = pd.concat([big, data[data[UNREG_COL] >= CALLOUT_MIN_CFS]]).drop_duplicates()

    for log_axes in (False, True):
        fig, ax = plt.subplots(figsize=(9.5, 8.5))
        ax.scatter(x, y, s=46, facecolor="#2c7fb8", edgecolor="0.25", lw=0.6,
                   zorder=3, label="Water year (adjusted regulated peak)")
        lim = [min(x.min(), y.min()) * 0.85, max(x.max(), y.max()) * 1.15]
        ax.plot(lim, lim, color="k", lw=1.2, ls="--", zorder=2,
                label="1:1 (no reservoir effect)")
        xs = np.linspace(max(lim[0], 1.0), lim[1], 200)
        ax.plot(xs, apply_power_law(fit, xs), color="#c0392b", lw=1.8, zorder=4,
                label="Fit: reg = %.4g x unreg$^{%.3f}$  (r$^2$=%.3f)"
                      % (fit["a"], fit["b"], fit["r2"]))
        band = 10 ** fit["se_dex"]
        ax.fill_between(xs, apply_power_law(fit, xs) / band,
                        apply_power_law(fit, xs) * band, color="#c0392b",
                        alpha=0.12, zorder=1, label="+/- 1 std error")
        annotate_points(ax, big[UNREG_COL].values, big["reg_peak"].values,
                        ["WY%d" % w for w in big["WY"]], log_axes)
        if log_axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.grid(which="both", alpha=0.25)
        else:
            ax.set_xlim(0, lim[1])
            ax.set_ylim(0, lim[1])
            ax.grid(alpha=0.3)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
        ax.set_xlabel("Unregulated peak at Castle Rock (cfs)")
        ax.set_ylabel("Adjusted regulated peak at Castle Rock (cfs)")
        ax.set_title("Castle Rock unregulated vs regulated peak%s\n"
                     "points below the 1:1 line are the reduction the project achieved"
                     % ("  (log-log)" if log_axes else ""), fontsize=11)
        ax.legend(loc="upper left", fontsize=9)
        fig.tight_layout()
        fig.savefig("%s_scatter_%s.png" % (stem, "loglog" if log_axes else "linear"),
                    dpi=150)
        plt.close(fig)


def probability_axis(ax, ticks, limits):
    """SSP-style normal-probability axis, AEP decreasing to the right."""
    z_ticks = stats.norm.ppf(1.0 - np.array(ticks))
    ax.set_xlim(stats.norm.ppf(1.0 - limits[0]), stats.norm.ppf(1.0 - limits[1]))
    ax.xaxis.set_major_locator(FixedLocator(z_ticks))
    ax.set_xticklabels(["%g" % (t * 100) if t * 100 >= 1 else "%g" % (t * 100)
                        for t in ticks])
    ax.set_xlabel("Annual exceedance probability (%)")


def plot_frequency(freq, data, fit, stem):
    """Unregulated and inferred regulated frequency curves, SSP idiom."""
    fig, ax = plt.subplots(figsize=(11, 8.5))

    z = stats.norm.ppf(1.0 - freq["AEP"].values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    reg_curve = apply_power_law(fit, unreg_curve)
    band = 10 ** fit["se_dex"]

    ax.plot(z, unreg_curve, color="#2c7fb8", lw=2.2, zorder=4,
            label="Unregulated peak (%s curve)" % FREQ_VALUE_COL.lower())
    ax.plot(z, reg_curve, color="#c0392b", lw=2.2, zorder=4,
            label="Regulated peak (inferred through the transform)")
    ax.fill_between(z, reg_curve / band, reg_curve * band, color="#c0392b",
                    alpha=0.14, zorder=1, label="Regulated, +/- 1 std error")

    if PLOTTING_BASIS == "from_curve":
        # AEP of each year's UNREGULATED peak, read off the unregulated curve;
        # the regulated peak of that same year is drawn at the same AEP.
        order = np.argsort(unreg_curve)
        aep_of = np.interp(data[UNREG_COL].values, unreg_curve[order],
                           freq["AEP"].values[order])
        zz = stats.norm.ppf(1.0 - np.clip(aep_of, 1e-6, 1 - 1e-6))
        ax.plot(zz, data[UNREG_COL].values, ls="none", marker="o", ms=5,
                mfc="#2c7fb8", mec="0.2", mew=0.5, zorder=5,
                label="Unregulated peaks (AEP from the curve)")
        ax.plot(zz, data["reg_peak"].values, ls="none", marker="s", ms=5,
                mfc="#c0392b", mec="0.2", mew=0.5, zorder=5,
                label="Adjusted regulated peaks (same AEP as their unreg peak)")
        for i in range(len(zz)):
            ax.plot([zz[i], zz[i]],
                    [data[UNREG_COL].values[i], data["reg_peak"].values[i]],
                    color="0.6", lw=0.5, zorder=3)
    else:
        uv, ua = weibull_plotting_positions(data[UNREG_COL].values)
        rv, ra = weibull_plotting_positions(data["reg_peak"].values)
        ax.plot(stats.norm.ppf(1.0 - ua), uv, ls="none", marker="o", ms=5,
                mfc="#2c7fb8", mec="0.2", mew=0.5, zorder=5,
                label="Unregulated peaks (Weibull, n=%d subset)" % len(uv))
        ax.plot(stats.norm.ppf(1.0 - ra), rv, ls="none", marker="s", ms=5,
                mfc="#c0392b", mec="0.2", mew=0.5, zorder=5,
                label="Adjusted regulated peaks (Weibull, n=%d subset)" % len(rv))

    # where the transform stops being supported by data
    supported = unreg_curve <= fit["x_max"]
    if supported.any() and not supported.all():
        z_edge = z[supported].max()
        ax.axvline(z_edge, color="0.35", lw=1.0, ls=":", zorder=2)
        ax.text(z_edge, FLOW_LIMITS[0] * 1.15,
                "  transform fitted to here\n  (unreg %s cfs)" % format(int(fit["x_max"]), ","),
                fontsize=8, color="0.3", va="bottom")

    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    probability_axis(ax, AEP_TICKS, AEP_LIMITS)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("Castle Rock peak flow frequency\n"
                 "regulated AEP inferred from the unregulated curve through "
                 "reg = %.4g x unreg$^{%.3f}$" % (fit["a"], fit["b"]), fontsize=12)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)
    fig.tight_layout()
    fig.savefig("%s_frequency.png" % stem, dpi=150)
    plt.close(fig)
    return reg_curve


def main():
    for path in (OUT_DIR, os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    data = pd.read_csv(DATASET_CSV).dropna(subset=[UNREG_COL, "reg_peak"])
    fit = fit_power_law(data[UNREG_COL].values, data["reg_peak"].values)

    freq = pd.read_csv(UNREG_FREQ_CSV)
    freq = freq[freq["Duration"] == FREQ_DURATION].sort_values("AEP", ascending=False)
    freq = freq.dropna(subset=["AEP", FREQ_VALUE_COL])

    print("=" * 78)
    print("Transform : reg = %.4g x unreg^%.4f   log r2 = %.4f   n = %d"
          % (fit["a"], fit["b"], fit["r2"], fit["n"]))
    print("Fitted over unregulated %s to %s cfs"
          % (format(int(fit["x_min"]), ","), format(int(fit["x_max"]), ",")))
    print("Scatter about the fit: x/ %.3f (1 sigma)" % 10 ** fit["se_dex"])
    print("=" * 78)

    plot_scatter(data, fit, PLOT_STEM)
    reg_curve = plot_frequency(freq, data, fit, PLOT_STEM)

    out = freq[["AEP", "Value", FREQ_VALUE_COL]].copy()
    out = out.rename(columns={"Value": "unreg_computed_cfs",
                              FREQ_VALUE_COL: "unreg_expected_cfs"})
    band = 10 ** fit["se_dex"]
    out["reg_inferred_cfs"] = reg_curve
    out["reg_lower_1se_cfs"] = reg_curve / band
    out["reg_upper_1se_cfs"] = reg_curve * band
    out["reduction_pct"] = 100.0 * (1.0 - reg_curve / out["unreg_expected_cfs"])
    out["extrapolated"] = out["unreg_expected_cfs"] > fit["x_max"]
    out.to_csv(os.path.join(OUT_DIR, "regulated_frequency_inferred.csv"),
               index=False, float_format="%.1f")

    show = out[out["AEP"].isin([0.5, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002])]
    print("Observed points placed by: %s" % PLOTTING_BASIS)
    print("\nINFERRED REGULATED CURVE")
    print(show[["AEP", "unreg_expected_cfs", "reg_inferred_cfs",
                "reg_lower_1se_cfs", "reg_upper_1se_cfs", "reduction_pct",
                "extrapolated"]].round(0).to_string(index=False))
    n_ex = int(out["extrapolated"].sum())
    if n_ex:
        print("\n   %d of %d AEPs are beyond the fitted range -- these depend on the"
              % (n_ex, len(out)))
        print("   power law holding past any event on record. The synthetics exist")
        print("   to replace that extrapolation with simulated points.")
    print("-" * 78)
    print("Plots : %s_scatter_linear.png, %s_scatter_loglog.png, %s_frequency.png"
          % (PLOT_STEM, PLOT_STEM, PLOT_STEM))
    print("Table : %s/regulated_frequency_inferred.csv" % OUT_DIR)


main()
