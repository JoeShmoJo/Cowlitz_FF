"""
Plot_2009_vs_2026_Frequency.py

Compare the 2009 Hydrology Restudy frequency curves with the 2026
CAS_Unreg_FF curves for the Peak, 1-, 3- and 5-day durations, on the
same HEC-SSP style Bulletin 17 axes used by
Plot_Duration_Frequency.py.

    2026 curves: solid
    2009 curves: dashed, same color per duration

DURATION MAPPING: the 2009 study's shortest duration is 0.20 day
(4.8 hr), not an instantaneous peak, so "Peak" pairs the 2026
instantaneous peak against the 2009 0.20-day statistics. Their adopted
means agree closely (4.784 vs 4.777), but the comparison is
approximate by construction. Columns 1, 3 and 5 map directly.

COMPUTED vs EXPECTED: the computed curve is the plain LP3 fit from the
adopted statistics. The expected probability curve adds Bulletin 17's
sampling-uncertainty adjustment, which depends on the RECORD LENGTH --
so it is not a property of the flood estimate but of how much data
stood behind it. The 2009 study had N = 69-79 and this one has
N = 93-95, so comparing expected curves would fold that ~20-year
difference into the comparison as if it were hydrology.

For an apples-to-apples comparison of the two studies, leave both
EP flags False: computed vs computed isolates the actual difference in
the fitted distributions. Set BOTH True only to reproduce whatever
each report published.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, pearson3, t as student_t

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
OUT_PNG = os.path.join(PROJECT_DIR, "output", "CAS_Unreg_2009_vs_2026.png")

TITLE = ("Cowlitz River at Castle Rock (14243000)\n"
         "2009 Hydrology Restudy vs 2026 Unregulated Frequency Curves")

# 2026 adopted LP3 statistics: label -> (log10 mean, log10 std, skew, N)
STATS_2026 = {
    "Peak":  (4.784, 0.202, -0.149, 95),
    "1-Day": (4.743, 0.203, -0.052, 94),
    "3-Day": (4.660, 0.192, -0.007, 93),
    "5-Day": (4.593, 0.179, -0.029, 93),
}

# 2009 Hydrology Restudy ADOPTED statistics, with systematic record
# length N. "Peak" uses the report's 0.20-day column (see note above).
STATS_2009 = {
    "Peak":  (4.777, 0.197, 0.193, 69),
    "1-Day": (4.742, 0.193, 0.180, 74),
    "3-Day": (4.663, 0.184, 0.168, 79),
    "5-Day": (4.597, 0.174, 0.155, 79),
}

COLORS = {"Peak": "#c00000", "1-Day": "#0070c0",
          "3-Day": "#00874e", "5-Day": "#7030a0"}

EP_2009 = False               # False = computed curve (apples-to-apples)
EP_2026 = False               # keep these two the SAME for comparison
SHOW_STATS_BOX = True

Y_LIMITS = (10_000, 300_000)       # flow axis, cfs
Y_TICKS = [10_000, 100_000, 300_000]
P_LIMITS = (0.999, 0.001)          # exceedance probability, left to right

P_TICKS = [0.999, 0.99, 0.9, 0.5, 0.2, 0.1,
           0.05, 0.02, 0.01, 0.005, 0.001]
T_TICKS = [1.001, 1.1, 2, 5, 10, 50, 200, 1000]
T_LABELS = ["1.0", "1.1", "2", "5", "10", "50", "200", "1000"]

CURVE_P_LIMITS = (0.99, 0.002)      # lines drawn 1.01-yr to 500-yr
CURVE_POINTS = 400
FIG_SIZE = (9.0, 7.2)
DPI = 200

###############################################################################
# FUNCTION DEFINITIONS


def prob_to_x(p):
    """Map exceedance probability to the normal-probability axis."""
    return norm.ppf(1.0 - np.asarray(p, dtype=float))


def curve_probabilities():
    """Probabilities spaced evenly on the probability axis."""
    x = np.linspace(prob_to_x(CURVE_P_LIMITS[0]),
                    prob_to_x(CURVE_P_LIMITS[1]), CURVE_POINTS)
    return 1.0 - norm.cdf(x)


def lp3_flows(mean_log, std_log, skew, probs):
    """Computed-curve flows for the given exceedance probabilities."""
    k = pearson3.ppf(1.0 - np.asarray(probs, dtype=float), skew=skew)
    return 10.0 ** (mean_log + k * std_log), k


def computed_prob_for_expected(p_expected, n):
    """Inverse of expected_probabilities: the computed exceedance
    probability whose EXPECTED probability is p_expected. Evaluating the
    LP3 curve here and plotting at p_expected draws the expected
    probability curve over exactly the requested probability range, so
    both curves span the same axis extent."""
    z = student_t.isf(np.asarray(p_expected, dtype=float), df=n - 1) \
        * np.sqrt(1.0 + 1.0 / n)
    return norm.sf(z)


def expected_probabilities(probs, n):
    """Bulletin 17 expected probability corresponding to each computed
    exceedance probability, for a record of length n. Normal theory:
        Pr[X > xbar + z*s] = Pr[t(n-1) > z / sqrt(1 + 1/n)]
    with z the NORMAL deviate of the computed probability. The
    adjustment is a function of probability and record length only --
    feeding the Pearson III K here instead would fold the skew into the
    shift and move the curve at the median, where the adjustment must
    vanish. The heavier t tail makes the expected probability exceed
    the computed one, so plotting each value at its expected
    probability lifts the curve in the tails, more so for short
    records."""
    z = norm.isf(np.asarray(probs, dtype=float))
    return student_t.sf(z / np.sqrt(1.0 + 1.0 / n), df=n - 1)


def duration_curve(stats, use_ep, probs):
    """(plot probabilities, flows) for one duration."""
    mean_log, std_log, skew, n = stats
    p_eval = computed_prob_for_expected(probs, n) if use_ep else probs
    q, _ = lp3_flows(mean_log, std_log, skew, p_eval)
    return probs, q


def draw_axes(ax):
    """Probability x axis, log flow y axis, return-period axis on top."""
    ax.set_xlim(prob_to_x(P_LIMITS[0]), prob_to_x(P_LIMITS[1]))
    ax.set_xticks(prob_to_x(P_TICKS))
    ax.set_xticklabels([("%g" % p) for p in P_TICKS], fontsize=8)
    ax.set_xlabel("Probability", fontsize=10)

    ax.set_yscale("log")
    ax.set_ylim(*Y_LIMITS)
    ax.set_yticks(Y_TICKS)
    ax.set_ylabel("Flow (cfs)", fontsize=10)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.tick_params(axis="y", labelsize=8)

    ax.grid(which="major", color="0.55", linewidth=0.6)
    ax.grid(which="minor", axis="y", color="0.80", linewidth=0.4)
    ax.set_axisbelow(True)

    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(prob_to_x([1.0 / v for v in T_TICKS]))
    top.set_xticklabels(T_LABELS, fontsize=8)
    top.set_xlabel("Return Period (years)", fontsize=10)
    return top


def stats_text():
    """Side-by-side summary of both statistic sets."""
    lines = [f"{'':<7}{'2009 adopted':>22}   {'2026 adopted':>22}",
             f"{'Dur':<7}{'Mean':>7}{'Std':>7}{'Skew':>8}   "
             f"{'Mean':>7}{'Std':>7}{'Skew':>8}"]
    for label in STATS_2026:
        a = STATS_2009[label]
        b = STATS_2026[label]
        lines.append(f"{label:<7}{a[0]:>7.3f}{a[1]:>7.3f}{a[2]:>8.3f}   "
                     f"{b[0]:>7.3f}{b[1]:>7.3f}{b[2]:>8.3f}")
    return "\n".join(lines)


def report(probs, curves_2009, curves_2026):
    """Console comparison at the standard probabilities."""
    targets = [0.5, 0.1, 0.02, 0.01, 0.002]
    print(f"{'Duration':<9}{'P':>8}{'2009':>11}{'2026':>11}"
          f"{'diff':>10}{'pct':>8}")
    for label in STATS_2026:
        p9, q9 = curves_2009[label]
        p6, q6 = curves_2026[label]
        for tp in targets:
            a = float(np.interp(tp, p9[::-1], q9[::-1]))
            b = float(np.interp(tp, p6[::-1], q6[::-1]))  # p9 == p6 == probs
            print(f"{label:<9}{tp:>8.3f}{a:>11,.0f}{b:>11,.0f}"
                  f"{b - a:>10,.0f}{100 * (b - a) / a:>7.1f}%")


###############################################################################
# MAIN


def main():
    probs = curve_probabilities()
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    draw_axes(ax)

    curves_2009, curves_2026 = {}, {}
    for label in STATS_2026:
        color = COLORS[label]

        p9, q9 = duration_curve(STATS_2009[label], EP_2009, probs)
        curves_2009[label] = (p9, q9)
        ax.plot(prob_to_x(p9), q9, color=color, linewidth=1.4,
                linestyle="--", dashes=(5, 3), zorder=3,
                label=f"{label} 2009"
                      + (" (expected)" if EP_2009 else ""))

        p6, q6 = duration_curve(STATS_2026[label], EP_2026, probs)
        curves_2026[label] = (p6, q6)
        ax.plot(prob_to_x(p6), q6, color=color, linewidth=1.7,
                linestyle="-", zorder=4,
                label=f"{label} 2026"
                      + (" (expected)" if EP_2026 else ""))

    report(probs, curves_2009, curves_2026)

    ax.set_title(TITLE, fontsize=11, pad=26)
    ax.legend(loc="upper left", fontsize=7.5, ncol=2, framealpha=0.9)

    if SHOW_STATS_BOX:
        ax.text(0.985, 0.03, stats_text(), transform=ax.transAxes,
                ha="right", va="bottom", family="monospace", fontsize=6.5,
                bbox=dict(boxstyle="square,pad=0.4", facecolor="white",
                          edgecolor="0.4", linewidth=0.6))

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=DPI)
    print(f"\nWrote {OUT_PNG}")


if __name__ == "__main__":
    main()
