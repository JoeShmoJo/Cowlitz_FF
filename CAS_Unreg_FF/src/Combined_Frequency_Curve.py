"""
Plot_Duration_Frequency.py

Plot the adopted Log-Pearson III frequency curves for all four durations
(Peak, 1-, 3-, 5-day) on one HEC-SSP style Bulletin 17 plot: log flow
axis, normal-probability x axis increasing in rarity to the right, and a
return-period axis across the top.

Curves are computed from the adopted LP3 statistics entered in STATS
below (log10 mean, log10 standard deviation, adopted skew), the same
values HEC-SSP reports:

    log10(Q) = mean + K(P, G) * std

K is the Pearson III frequency factor, taken from the exact Pearson III
quantile function rather than the Wilson-Hilferty approximation.

Observed events are optionally overlaid from the assembled record
(wy_record_ssp.csv) using Weibull plotting positions. NOTE: HEC-SSP
plots Hirsch-Stedinger positions, so the points here will sit slightly
differently on the probability axis; the fitted curves are unaffected.
Set PLOT_OBSERVED = False for curves only.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, pearson3

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
RECORD_CSV = os.path.join(PROJECT_DIR, "output", "wy_record_ssp.csv")
OUT_PNG = os.path.join(PROJECT_DIR, "output", "CAS_Unreg_duration_frequency.png")

TITLE = "Cowlitz River at Castle Rock (14243000) — Unregulated Frequency Curves"

# Adopted LP3 statistics (log10 units), one entry per duration:
#   label -> (mean, std, adopted skew, record column, color)
# The record column is the matching header in wy_record_ssp.csv and is
# only used when PLOT_OBSERVED is True.
STATS = {
    "Peak":  (4.784, 0.202, -0.149, "Peak",      "#c00000"),
    "1-Day": (4.743, 0.203, -0.052, "One_day",   "#0070c0"),
    "3-Day": (4.660, 0.192, -0.007, "Three_Day", "#00874e"),
    "5-Day": (4.593, 0.179, -0.029, "Five_Day",  "#7030a0"),
}

PLOT_OBSERVED = True          # overlay observed events as open circles
SHOW_STATS_BOX = True         # SSP-style statistics box

Y_LIMITS = (10_000, 1_000_000)     # flow axis, cfs
P_LIMITS = (0.9999, 0.0001)        # exceedance probability, left to right

# Bottom axis (exceedance probability) and top axis (return period, yrs)
P_TICKS = [0.9999, 0.999, 0.99, 0.9, 0.5, 0.2, 0.1,
           0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0001]
T_TICKS = [1.0001, 1.1, 2, 5, 10, 50, 200, 1000, 10000]
T_LABELS = ["1.0", "1.1", "2", "5", "10", "50", "200", "1000", "10000"]

# Fitted curves are independent LP3 fits, so nothing forces the shorter
# duration to stay above the longer one at every probability. Where the
# adopted skews differ enough, the curves can cross in the far tail
# (Peak vs 1-Day cross near T = 5,400 yr with the current statistics).
# The check below reports any crossing; it does not alter the curves.
CHECK_CURVE_ORDER = True

CURVE_POINTS = 400            # resolution of the plotted curves
FIG_SIZE = (9.0, 7.0)
DPI = 200

###############################################################################
# FUNCTION DEFINITIONS


def prob_to_x(p):
    """Map exceedance probability to the normal-probability axis.
    Rare events (small p) plot to the right, as in HEC-SSP."""
    return norm.ppf(1.0 - np.asarray(p, dtype=float))


def lp3_curve(mean_log, std_log, skew, probs):
    """Log-Pearson III flows for the given exceedance probabilities.
    K is the exact Pearson III frequency factor for the adopted skew."""
    k = pearson3.ppf(1.0 - np.asarray(probs, dtype=float), skew=skew)
    return 10.0 ** (mean_log + k * std_log)


def curve_probabilities():
    """Probabilities spaced evenly on the probability axis, so curves are
    smooth across the whole plot rather than bunched near the median."""
    x = np.linspace(prob_to_x(P_LIMITS[0]), prob_to_x(P_LIMITS[1]),
                    CURVE_POINTS)
    return 1.0 - norm.cdf(x)


def observed_points(column):
    """Observed values and Weibull plotting positions for one duration.
    Returns (probabilities, sorted flows); empty arrays if unavailable."""
    empty = (np.array([]), np.array([]))
    if not os.path.exists(RECORD_CSV):
        return empty
    import pandas as pd
    df = pd.read_csv(RECORD_CSV)
    if column not in df.columns:
        print(f"  column {column} not in {os.path.basename(RECORD_CSV)}")
        return empty
    v = df[column].to_numpy(dtype=float)
    v = np.sort(v[np.isfinite(v) & (v > 0)])[::-1]   # largest first
    if v.size == 0:
        return empty
    n = v.size
    ranks = np.arange(1, n + 1)
    return ranks / (n + 1.0), v          # Weibull: p = i / (n + 1)


def draw_axes(ax):
    """Probability x axis, log flow y axis, return-period axis on top."""
    ax.set_xlim(prob_to_x(P_LIMITS[0]), prob_to_x(P_LIMITS[1]))
    ax.set_xticks(prob_to_x(P_TICKS))
    ax.set_xticklabels([("%g" % p) for p in P_TICKS], fontsize=8)
    ax.set_xlabel("Probability", fontsize=10)

    ax.set_yscale("log")
    ax.set_ylim(*Y_LIMITS)
    ax.set_ylabel("Flow (cfs)", fontsize=10)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.tick_params(axis="y", labelsize=8)

    ax.grid(which="major", color="0.55", linewidth=0.6)
    ax.grid(which="minor", axis="y", color="0.80", linewidth=0.4)
    ax.set_axisbelow(True)

    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(prob_to_x([1.0 / t for t in T_TICKS]))
    top.set_xticklabels(T_LABELS, fontsize=8)
    top.set_xlabel("Return Period (years)", fontsize=10)
    return top


def check_curve_order(curves):
    """Report probabilities where a shorter duration falls below a longer
    one. Physically impossible, but independent LP3 fits do not prevent
    it -- differing adopted skews can cross the curves in the far tail."""
    labels = list(curves.keys())
    problems = []
    for short, long_ in zip(labels[:-1], labels[1:]):
        qs, ql = curves[short][1], curves[long_][1]
        bad = qs < ql
        if bad.any():
            p = curves[short][0][bad]
            problems.append((short, long_, p.max(), p.min(),
                             float(np.max(ql[bad] - qs[bad]))))
    if not problems:
        print("\nCurve order OK: each duration stays above the next-longer "
              "one across the plotted range.")
        return
    print("\n*** CURVE CROSSING (shorter duration falls below longer) ***")
    for short, long_, p_hi, p_lo, gap in problems:
        print(f"  {short} < {long_} for P <= {p_hi:.5f} "
              f"(T >= {1.0 / p_hi:,.0f} yr), max gap {gap:,.0f} cfs")
    print("  Independent LP3 fits do not enforce ordering between "
          "durations; review the adopted skews if this falls inside the "
          "range you intend to report.")


def stats_text():
    """SSP-style summary of the adopted statistics."""
    lines = [f"{'Duration':<8}{'Mean':>7}{'Std':>7}{'Skew':>8}"]
    for label, (m, s, g, _, _) in STATS.items():
        lines.append(f"{label:<8}{m:>7.3f}{s:>7.3f}{g:>8.3f}")
    return "\n".join(lines)


###############################################################################
# MAIN


def main():
    probs = curve_probabilities()
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    draw_axes(ax)

    curves = {}
    for label, (mean_log, std_log, skew, column, color) in STATS.items():
        q = lp3_curve(mean_log, std_log, skew, probs)
        curves[label] = (probs, q)
        ax.plot(prob_to_x(probs), q, color=color, linewidth=1.6,
                label=f"{label} computed curve", zorder=4)
        print(f"{label:<6} 50%={q[np.argmin(np.abs(probs - 0.5))]:>9,.0f}  "
              f"10%={q[np.argmin(np.abs(probs - 0.1))]:>9,.0f}  "
              f"1%={q[np.argmin(np.abs(probs - 0.01))]:>9,.0f}  "
              f"0.2%={q[np.argmin(np.abs(probs - 0.002))]:>9,.0f} cfs")

        if PLOT_OBSERVED:
            p_obs, q_obs = observed_points(column)
            if q_obs.size:
                ax.plot(prob_to_x(p_obs), q_obs, linestyle="none",
                        marker="o", markersize=3.2, markerfacecolor="none",
                        markeredgecolor=color, markeredgewidth=0.7,
                        label=f"{label} observed (n={q_obs.size})", zorder=3)

    if CHECK_CURVE_ORDER:
        check_curve_order(curves)

    ax.set_title(TITLE, fontsize=11, pad=28)
    ax.legend(loc="upper left", fontsize=7.5, ncol=2, framealpha=0.9)

    if SHOW_STATS_BOX:
        ax.text(0.985, 0.03, stats_text(), transform=ax.transAxes,
                ha="right", va="bottom", family="monospace", fontsize=7,
                bbox=dict(boxstyle="square,pad=0.4", facecolor="white",
                          edgecolor="0.4", linewidth=0.6))

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=DPI)
    print(f"\nWrote {OUT_PNG}")


if __name__ == "__main__":
    main()