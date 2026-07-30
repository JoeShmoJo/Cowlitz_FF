"""
Frequency_Curves_And_Table.py

One script, three outputs:

  1. output/CAS_Unreg_2026_frequency.png   2026 curves only
  2. output/CAS_Unreg_2009_vs_2026.png     2026 vs 2009 Hydrology Restudy
  3. output/CAS_Unreg_frequency_table.csv  appendix table, 2026 results

PLOTS are drawn as Log-Pearson III curves from the adopted statistics
(log10 mean, std, adopted skew) on HEC-SSP style Bulletin 17 axes. That
reproduces SSP's Computed Curve to within 0.32% across the full range;
the small residual is only because the adopted statistics are reported
to three decimals. Computing both studies the same way keeps the
comparison apples-to-apples -- see EXPECTED PROBABILITY below.

The APPENDIX TABLE is parsed straight out of the SSP report files
rather than recomputed, because Variance Log, the expected probability
curve and the EMA confidence limits are products of the EMA fit and
cannot be recovered from the three adopted moments.

EXPECTED PROBABILITY: SSP reports a Computed Curve and an Expected
Probability curve. The latter adds Bulletin 17's sampling-uncertainty
adjustment, which is a function of probability and RECORD LENGTH, not
of the flood estimate. The 2009 study had N = 69-79 against 93-95 here,
so comparing expected curves would fold that record-length difference
into the comparison as if it were hydrology. Both plots therefore use
computed curves. The appendix table carries SSP's expected probability
column alongside the computed one, so both are on the record.

DURATION MAPPING for the comparison: the 2009 study's shortest duration
is 0.20 day (4.8 hr), not an instantaneous peak, so "Peak" pairs the
2026 instantaneous peak against the 2009 0.20-day statistics. Their
adopted means agree closely (4.784 vs 4.777), but the pairing is
approximate by construction. Columns 1, 3 and 5 map directly.
"""

import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, pearson3

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
OUT_PNG_2026 = os.path.join(PROJECT_DIR, "output",
                            "CAS_Unreg_2026_frequency.png")
OUT_PNG_COMPARE = os.path.join(PROJECT_DIR, "output",
                               "CAS_Unreg_2009_vs_2026.png")
OUT_TABLE_CSV = os.path.join(PROJECT_DIR, "output",
                             "CAS_Unreg_frequency_table.csv")

# SSP Bulletin 17 report files, one per duration. Check the analysis
# names against your Bulletin17Results folder.
SSP_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF", "ssp", "2026_Restudy", "Bulletin17Results")
SSP_REPORTS = {
    "Peak":  os.path.join(SSP_DIR, "CAS_2026_p", "CAS_2026_p.rpt"),
    "1-Day": os.path.join(SSP_DIR, "CAS_2026_1", "CAS_2026_1.rpt"),
    "3-Day": os.path.join(SSP_DIR, "CAS_2026_3", "CAS_2026_3.rpt"),
    "5-Day": os.path.join(SSP_DIR, "CAS_2026_5", "CAS_2026_5.rpt"),
}

LOCATION_NAME = "Castle Rock"
DURATION_DAYS = {"Peak": 0, "1-Day": 1, "3-Day": 3, "5-Day": 5}

# 2026 adopted LP3 statistics: label -> (log10 mean, log10 std, skew)
STATS_2026 = {
    "Peak":  (4.778, 0.202, -0.135),
    "1-Day": (4.743, 0.203, -0.052),
    "3-Day": (4.660, 0.192, -0.007),
    "5-Day": (4.593, 0.179, -0.029),
}

# 2009 Hydrology Restudy ADOPTED statistics. "Peak" is the report's
# 0.20-day column (see note above).
STATS_2009 = {
    "Peak":  (4.777, 0.197, 0.193),
    "1-Day": (4.742, 0.193, 0.180),
    "3-Day": (4.663, 0.184, 0.168),
    "5-Day": (4.597, 0.174, 0.155),
}

COLORS = {"Peak": "#c00000", "1-Day": "#0070c0",
          "3-Day": "#00874e", "5-Day": "#7030a0"}

SHOW_STATS_BOX = True

Y_LIMITS = (10_000, 300_000)       # flow axis, cfs
Y_TICKS = [10_000, 100_000, 300_000]
P_LIMITS = (0.999, 0.001)          # axis extent, exceedance probability

P_TICKS = [0.999, 0.99, 0.9, 0.5, 0.2, 0.1,
           0.05, 0.02, 0.01, 0.005, 0.001]
T_TICKS = [1.001, 1.1, 2, 5, 10, 50, 200, 1000]
T_LABELS = ["1.0", "1.1", "2", "5", "10", "50", "200", "1000"]

CURVE_P_LIMITS = (0.99, 0.002)      # lines drawn 1.01-yr to 500-yr
CURVE_POINTS = 400
FIG_SIZE = (9.0, 7.0)
DPI = 200

###############################################################################
# FUNCTION DEFINITIONS


def prob_to_x(p):
    """Map exceedance probability to the normal-probability axis."""
    return norm.ppf(1.0 - np.asarray(p, dtype=float))


def curve_probabilities():
    """Probabilities spaced evenly along the probability axis."""
    x = np.linspace(prob_to_x(CURVE_P_LIMITS[0]),
                    prob_to_x(CURVE_P_LIMITS[1]), CURVE_POINTS)
    return 1.0 - norm.cdf(x)


def lp3_flows(mean_log, std_log, skew, probs):
    """Computed-curve flows for the given exceedance probabilities."""
    k = pearson3.ppf(1.0 - np.asarray(probs, dtype=float), skew=skew)
    return 10.0 ** (mean_log + k * std_log)


def parse_ssp_report(path):
    """Pull the Frequency Curve table out of an SSP .rpt file.

    Data rows look like:
      |  126,877.9   0.00367377   136,163.4 |   0.2  |  166,906.7  105,871.8 |
       computed     variance log  expected     pct       0.05 limit  0.95 limit

    Returns a DataFrame [pct_chance, computed, variance_log, expected,
    conf_05, conf_95], or an empty frame if the file is unreadable.
    """
    cols = ["pct_chance", "computed", "variance_log", "expected",
            "conf_05", "conf_95"]
    if not os.path.exists(path):
        print(f"  MISSING report: {path}")
        return pd.DataFrame(columns=cols)

    num = re.compile(r"-?[\d,]+\.?\d*(?:[eE][+-]?\d+)?")
    rows, in_block = [], False
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if "Frequency Curve" in line:
                in_block = True
                continue
            if not in_block:
                continue
            parts = line.split("|")
            if len(parts) < 4:
                if rows:            # table finished
                    break
                continue
            left = num.findall(parts[1])
            mid = num.findall(parts[2])
            right = num.findall(parts[3])
            if len(left) != 3 or len(mid) != 1 or len(right) != 2:
                continue
            f = lambda s: float(s.replace(",", ""))
            rows.append({
                "pct_chance": f(mid[0]),
                "computed": f(left[0]),
                "variance_log": f(left[1]),
                "expected": f(left[2]),
                "conf_05": f(right[0]),
                "conf_95": f(right[1]),
            })
    if not rows:
        print(f"  no Frequency Curve rows found in {path}")
    return pd.DataFrame(rows, columns=cols)


def build_appendix_table():
    """Assemble the appendix CSV from the SSP reports, in the agreed
    column layout. AEP is a fraction; UpperConf is SSP's 0.05 limit and
    LowerConf its 0.95 limit."""
    out = []
    for label, path in SSP_REPORTS.items():
        print(f"Reading {label}: {path}")
        df = parse_ssp_report(path)
        if df.empty:
            continue
        print(f"  {len(df)} ordinates")
        aep = df["pct_chance"] / 100.0
        out.append(pd.DataFrame({
            "LocationName": LOCATION_NAME,
            "Duration": label,
            "DurationDays": DURATION_DAYS.get(label, ""),
            "AEP": aep,
            "VarianceLog": df["variance_log"],
            "zScore": norm.ppf(1.0 - aep),
            "Value": df["computed"],
            "UpperConf": df["conf_05"],
            "LowerConf": df["conf_95"],
            "Expected": df["expected"],
        }))
    if not out:
        print("No SSP reports parsed; appendix table not written.")
        return pd.DataFrame()
    table = pd.concat(out, ignore_index=True).sort_values(
        ["DurationDays", "AEP"], kind="stable")
    os.makedirs(os.path.dirname(OUT_TABLE_CSV), exist_ok=True)
    table.to_csv(OUT_TABLE_CSV, index=False)
    print(f"\nWrote {len(table)} rows -> {OUT_TABLE_CSV}")
    return table


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


def stats_box_text(include_2009):
    """Monospace statistics summary for the plot corner."""
    if not include_2009:
        lines = [f"{'Duration':<9}{'Mean':>7}{'Std':>7}{'Skew':>8}"]
        for label, (m, s, g) in STATS_2026.items():
            lines.append(f"{label:<9}{m:>7.3f}{s:>7.3f}{g:>8.3f}")
        return "\n".join(lines)
    lines = [f"{'':<7}{'2009 adopted':>22}   {'2026 adopted':>22}",
             f"{'Dur':<7}{'Mean':>7}{'Std':>7}{'Skew':>8}   "
             f"{'Mean':>7}{'Std':>7}{'Skew':>8}"]
    for label in STATS_2026:
        a, b = STATS_2009[label], STATS_2026[label]
        lines.append(f"{label:<7}{a[0]:>7.3f}{a[1]:>7.3f}{a[2]:>8.3f}   "
                     f"{b[0]:>7.3f}{b[1]:>7.3f}{b[2]:>8.3f}")
    return "\n".join(lines)


def make_plot(include_2009, out_png, title):
    """Draw the frequency curves; optionally overlay the 2009 study."""
    probs = curve_probabilities()
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    draw_axes(ax)

    for label, (mean_log, std_log, skew) in STATS_2026.items():
        color = COLORS[label]
        if include_2009:
            a = STATS_2009[label]
            ax.plot(prob_to_x(probs), lp3_flows(*a, probs), color=color,
                    linewidth=1.4, linestyle="--", dashes=(5, 3), zorder=3,
                    label=f"{label} 2009")
        ax.plot(prob_to_x(probs), lp3_flows(mean_log, std_log, skew, probs),
                color=color, linewidth=1.7, zorder=4,
                label=f"{label} 2026" if include_2009 else label)

    ax.set_title(title, fontsize=11, pad=26)
    ax.legend(loc="upper left", fontsize=7.5,
              ncol=2 if include_2009 else 1, framealpha=0.9)

    if SHOW_STATS_BOX:
        ax.text(0.985, 0.03, stats_box_text(include_2009),
                transform=ax.transAxes, ha="right", va="bottom",
                family="monospace", fontsize=6.5 if include_2009 else 7,
                bbox=dict(boxstyle="square,pad=0.4", facecolor="white",
                          edgecolor="0.4", linewidth=0.6))

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Wrote {out_png}")


def report_comparison():
    """Console summary of 2026 vs 2009 at the standard probabilities."""
    targets = [0.5, 0.1, 0.02, 0.01, 0.002]
    print(f"\n{'Duration':<9}{'AEP':>8}{'2009':>11}{'2026':>11}"
          f"{'diff':>10}{'pct':>8}")
    for label in STATS_2026:
        a = lp3_flows(*STATS_2009[label], np.array(targets))
        b = lp3_flows(*STATS_2026[label], np.array(targets))
        for tp, x, y in zip(targets, a, b):
            print(f"{label:<9}{tp:>8.3f}{x:>11,.0f}{y:>11,.0f}"
                  f"{y - x:>10,.0f}{100 * (y - x) / x:>7.1f}%")


###############################################################################
# MAIN


def main():
    make_plot(False, OUT_PNG_2026,
              "Cowlitz River at Castle Rock (14243000)\n"
              "2026 Unregulated Frequency Curves")
    make_plot(True, OUT_PNG_COMPARE,
              "Cowlitz River at Castle Rock (14243000)\n"
              "2009 Hydrology Restudy vs 2026 Unregulated Frequency Curves")
    report_comparison()
    print()
    build_appendix_table()


if __name__ == "__main__":
    main()
