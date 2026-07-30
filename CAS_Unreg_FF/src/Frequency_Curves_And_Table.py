"""
Frequency_Curves_And_Table.py

Reads the HEC-SSP Bulletin 17 results committed under CAS_Unreg_FF/ssp
and produces every frequency deliverable for the memo:

  Section 10 figures
    output/CAS_Unreg_2026_frequency.png    2026 curves, all four durations
    output/CAS_Unreg_2009_vs_2026.png      the same with the 2009 study

  Section 10 table
    output/CAS_Unreg_distribution_parameters.csv
        distribution parameters and event counts per duration

  Appendix figures (one per duration, SSP layout)
    output/CAS_Unreg_freq_Peak.png  ... _1-Day, _3-Day, _5-Day
        computed curve, expected probability curve, 5%/95% confidence
        limits, observed events at Hirsch-Stedinger plotting positions

  Appendix table
    output/CAS_Unreg_frequency_table.csv

ONLY the 2009 statistics are hard coded. Every 2026 number -- moments,
skews, event counts, curve ordinates, confidence limits and plotting
positions -- is parsed from the SSP reports, so re-running SSP and then
this script keeps the memo in step.

ADOPTED ANALYSES: the peak uses CAS_2026_p_Sensor_1969_1973, in which
WY1969-1973 are censored (entered as intervals rather than values) and
no low-outlier test is applied. The duration analyses need no censoring
because those water years carry no duration data at all.

COMPUTED vs EXPECTED: SSP reports both. The expected probability curve
adds Bulletin 17's sampling-uncertainty adjustment, a function of
probability and record length rather than of the flood estimate. The
2009 study had shorter records, so the study comparison uses computed
curves on both sides; the per-duration appendix figures show both, as
SSP does.

DURATION MAPPING for the comparison: the 2009 study's shortest duration
is 0.20 day (4.8 hr), not an instantaneous peak, so "Peak" pairs the
2026 instantaneous peak against the 2009 0.20-day statistics. Their
means agree closely (4.778 vs 4.777), but the pairing is approximate.
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
OUT_DIR = os.path.join(PROJECT_DIR, "output")
SSP_DIR = os.path.join(PROJECT_DIR, "ssp", "2026_Restudy", "Bulletin17Results")

# Adopted SSP analysis per duration (folder name = report stem).
SSP_ANALYSES = {
    "Peak":  "CAS_2026_p_Sensor_1969_1973",   # WY1969-1973 censored, no PILF
    "1-Day": "CAS_2026_1",
    "3-Day": "CAS_2026_3",
    "5-Day": "CAS_2026_5",
}

# 2009 Hydrology Restudy ADOPTED statistics (hard coded; the only
# non-SSP numbers here). "Peak" is the report's 0.20-day column.
STATS_2009 = {
    "Peak":  (4.777, 0.197, 0.193),      # log10 mean, log10 std, skew
    "1-Day": (4.742, 0.193, 0.180),
    "3-Day": (4.663, 0.184, 0.168),
    "5-Day": (4.597, 0.174, 0.155),
}

LOCATION_NAME = "Castle Rock"
DURATION_DAYS = {"Peak": 0, "1-Day": 1, "3-Day": 3, "5-Day": 5}
COLORS = {"Peak": "#c00000", "1-Day": "#0070c0",
          "3-Day": "#00874e", "5-Day": "#7030a0"}

SHOW_STATS_BOX = True

# Summary figures (Section 10)
Y_LIMITS = (10_000, 300_000)
Y_TICKS = [10_000, 100_000, 300_000]
CURVE_P_LIMITS = (0.99, 0.002)      # lines drawn 1.01-yr to 500-yr

# Per-duration figures (Appendix) -- taller, to hold the confidence band
Y_LIMITS_DUR = (10_000, 1_000_000)

P_LIMITS = (0.999, 0.001)           # axis extent, exceedance probability
P_TICKS = [0.999, 0.99, 0.9, 0.5, 0.2, 0.1,
           0.05, 0.02, 0.01, 0.005, 0.001]
T_TICKS = [1.001, 1.1, 2, 5, 10, 50, 200, 1000]
T_LABELS = ["1.0", "1.1", "2", "5", "10", "50", "200", "1000"]

CURVE_POINTS = 400
FIG_SIZE = (9.0, 7.0)
DPI = 200

PARAM_ROWS = ["Mean", "Standard Dev", "Station Skew", "Regional Skew",
              "Weighted Skew", "Adopted Skew"]
DIAG_ROWS = ["EMA MSE [G at-site]", "Grubbs-Beck Critical Value",
             "Equivalent Record Length"]
EVENT_ROWS = ["Systematic Events", "Historic Events", "High Outliers",
              "Low Outliers", "Zero Events", "Missing Events",
              "Historic Period"]

###############################################################################
# SSP REPORT PARSING


def report_path(label):
    stem = SSP_ANALYSES[label]
    return os.path.join(SSP_DIR, stem, stem + ".rpt")


def read_report(label):
    path = report_path(label)
    if not os.path.exists(path):
        print(f"  MISSING report: {path}")
        return ""
    with open(path, "r", errors="ignore") as fh:
        return fh.read()


def parse_frequency_curve(txt):
    """Frequency Curve block -> DataFrame."""
    cols = ["pct_chance", "computed", "variance_log", "expected",
            "conf_05", "conf_95"]
    if "<< Frequency Curve >>" not in txt:
        return pd.DataFrame(columns=cols)
    blk = txt.split("<< Frequency Curve >>")[1].split("<<")[0]
    num = re.compile(r"-?[\d,]+\.?\d*(?:[eE][+-]?\d+)?")
    rows = []
    for line in blk.splitlines():
        parts = line.split("|")
        if len(parts) < 4:
            continue
        left, mid, right = (num.findall(parts[1]), num.findall(parts[2]),
                            num.findall(parts[3]))
        if len(left) != 3 or len(mid) != 1 or len(right) != 2:
            continue
        f = lambda s: float(s.replace(",", ""))
        rows.append({"pct_chance": f(mid[0]), "computed": f(left[0]),
                     "variance_log": f(left[1]), "expected": f(left[2]),
                     "conf_05": f(right[0]), "conf_95": f(right[1])})
    return pd.DataFrame(rows, columns=cols)


def parse_statistics(txt):
    """Systematic Statistics block plus fitted-moment diagnostics."""
    out = {}
    if "<< Systematic Statistics >>" in txt:
        blk = txt.split("<< Systematic Statistics >>")[1]
        for line in blk.splitlines():
            if "|" not in line:
                continue
            for half in [p.strip() for p in line.split("|") if p.strip()]:
                m = re.match(r"^([A-Za-z][A-Za-z0-9 \-\[\]\(\)\.']*?)"
                             r"\s+(-?[\d\.]+)$", half)
                if m:
                    out[m.group(1).strip()] = float(m.group(2))
    for key, pat in [
            ("EMA MSE [G at-site]",
             r"EMA Estimate of MSE\[G at-site\]\s+([\d\.\-]+)"),
            ("Equivalent Record Length",
             r"Equivalent Record Length \[Cohn et al \(1997\)\]"
             r"\s+([\d\.\-]+)"),
            ("Grubbs-Beck Critical Value",
             r"Grubbs-Beck Critical Value\s+([\d\.\-]+)")]:
        m = re.search(pat, txt)
        if m:
            out[key] = float(m.group(1))
    m = re.search(r"EMA w/ regional info and B17b MSE\(G\)\s+"
                  r"([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)",
                  txt)
    if m:
        out["_mean"] = float(m.group(1))
        out["_std"] = float(m.group(3))
        out["_skew"] = float(m.group(4))
    return out


def parse_plotting_positions(txt):
    """Ordered Events -> DataFrame [rank, wy, flow, pct_chance]."""
    cols = ["rank", "wy", "flow", "pct_chance"]
    if "<< Plotting Positions >>" not in txt:
        return pd.DataFrame(columns=cols)
    blk = txt.split("<< Plotting Positions >>")[1].split("<<")[0]
    rows = []
    for line in blk.splitlines():
        parts = line.split("|")
        if len(parts) < 3:
            continue
        m = re.findall(r"[\d,]+\.?\d*", parts[2])
        if len(m) == 4:
            rows.append({"rank": int(m[0]), "wy": int(m[1]),
                         "flow": float(m[2].replace(",", "")),
                         "pct_chance": float(m[3])})
    return pd.DataFrame(rows, columns=cols)


def load_all():
    data = {}
    for label in SSP_ANALYSES:
        txt = read_report(label)
        if not txt:
            continue
        data[label] = {"curve": parse_frequency_curve(txt),
                       "stats": parse_statistics(txt),
                       "plotpos": parse_plotting_positions(txt)}
        s = data[label]["stats"]
        print(f"{label:<6} mean {s.get('_mean', np.nan):.6f}  "
              f"std {s.get('_std', np.nan):.6f}  "
              f"skew {s.get('_skew', np.nan):+.6f}  "
              f"n={int(s.get('Systematic Events', 0))}  "
              f"{len(data[label]['curve'])} ordinates  "
              f"{len(data[label]['plotpos'])} events")
    return data


###############################################################################
# CURVES AND AXES


def prob_to_x(p):
    return norm.ppf(1.0 - np.asarray(p, dtype=float))


def curve_probabilities(limits):
    x = np.linspace(prob_to_x(limits[0]), prob_to_x(limits[1]), CURVE_POINTS)
    return 1.0 - norm.cdf(x)


def lp3_flows(mean_log, std_log, skew, probs):
    k = pearson3.ppf(1.0 - np.asarray(probs, dtype=float), skew=skew)
    return 10.0 ** (mean_log + k * std_log)


def draw_axes(ax, y_limits, y_ticks=None):
    ax.set_xlim(prob_to_x(P_LIMITS[0]), prob_to_x(P_LIMITS[1]))
    ax.set_xticks(prob_to_x(P_TICKS))
    ax.set_xticklabels([("%g" % p) for p in P_TICKS], fontsize=8)
    ax.set_xlabel("Probability", fontsize=10)

    ax.set_yscale("log")
    ax.set_ylim(*y_limits)
    if y_ticks:
        ax.set_yticks(y_ticks)
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


###############################################################################
# FIGURES


def stats_box_text(data, include_2009):
    if not include_2009:
        lines = [f"{'Duration':<9}{'Mean':>7}{'Std':>7}{'Skew':>8}{'n':>5}"]
        for label in SSP_ANALYSES:
            if label not in data:
                continue
            s = data[label]["stats"]
            lines.append(f"{label:<9}{s['_mean']:>7.3f}{s['_std']:>7.3f}"
                         f"{s['_skew']:>8.3f}"
                         f"{int(s.get('Systematic Events', 0)):>5d}")
        return "\n".join(lines)
    lines = [f"{'':<7}{'2009 adopted':>22}   {'2026 adopted':>22}",
             f"{'Dur':<7}{'Mean':>7}{'Std':>7}{'Skew':>8}   "
             f"{'Mean':>7}{'Std':>7}{'Skew':>8}"]
    for label in SSP_ANALYSES:
        if label not in data:
            continue
        a = STATS_2009[label]
        s = data[label]["stats"]
        lines.append(f"{label:<7}{a[0]:>7.3f}{a[1]:>7.3f}{a[2]:>8.3f}   "
                     f"{s['_mean']:>7.3f}{s['_std']:>7.3f}{s['_skew']:>8.3f}")
    return "\n".join(lines)


def make_summary_plot(data, include_2009, out_png, title):
    probs = curve_probabilities(CURVE_P_LIMITS)
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    draw_axes(ax, Y_LIMITS, Y_TICKS)

    for label in SSP_ANALYSES:
        if label not in data:
            continue
        color = COLORS[label]
        if include_2009:
            ax.plot(prob_to_x(probs), lp3_flows(*STATS_2009[label], probs),
                    color=color, linewidth=1.4, linestyle="--",
                    dashes=(5, 3), zorder=3, label=f"{label} 2009")
        s = data[label]["stats"]
        ax.plot(prob_to_x(probs),
                lp3_flows(s["_mean"], s["_std"], s["_skew"], probs),
                color=color, linewidth=1.7, zorder=4,
                label=f"{label} 2026" if include_2009 else label)

    ax.set_title(title, fontsize=11, pad=26)
    ax.legend(loc="upper left", fontsize=7.5,
              ncol=2 if include_2009 else 1, framealpha=0.9)
    if SHOW_STATS_BOX:
        ax.text(0.985, 0.03, stats_box_text(data, include_2009),
                transform=ax.transAxes, ha="right", va="bottom",
                family="monospace", fontsize=6.5 if include_2009 else 7,
                bbox=dict(boxstyle="square,pad=0.4", facecolor="white",
                          edgecolor="0.4", linewidth=0.6))
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Wrote {out_png}")


def make_duration_plot(data, label, out_png):
    """SSP-style single-duration figure."""
    d = data[label]
    curve, stats, pp = d["curve"], d["stats"], d["plotpos"]
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    draw_axes(ax, Y_LIMITS_DUR)

    probs = curve_probabilities(P_LIMITS)
    ax.plot(prob_to_x(probs),
            lp3_flows(stats["_mean"], stats["_std"], stats["_skew"], probs),
            color="#c00000", linewidth=1.6, zorder=5,
            label="Computed curve")

    if not curve.empty:
        p = (curve["pct_chance"] / 100.0).values
        order = np.argsort(-p)
        px = prob_to_x(p[order])
        ax.plot(px, curve["expected"].values[order], color="#0000c0",
                linewidth=1.2, linestyle=":", zorder=4,
                label="Expected probability curve")
        ax.plot(px, curve["conf_05"].values[order], color="#008000",
                linewidth=1.0, linestyle="--", dashes=(6, 3), zorder=3,
                label="5% confidence limit")
        ax.plot(px, curve["conf_95"].values[order], color="#008000",
                linewidth=1.0, linestyle="--", dashes=(2, 2), zorder=3,
                label="95% confidence limit")

    if not pp.empty:
        ax.plot(prob_to_x(pp["pct_chance"] / 100.0), pp["flow"],
                linestyle="none", marker="o", markersize=3.6,
                markerfacecolor="none", markeredgecolor="#0000c0",
                markeredgewidth=0.8, zorder=6,
                label="Observed events (Hirsch-Stedinger)")

    ax.set_title("Cowlitz River at Castle Rock (14243000)\n"
                 f"Unregulated {label} Frequency Curve — "
                 f"{SSP_ANALYSES[label]}", fontsize=11, pad=26)
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9)

    box = (f"{'Mean':<21}{stats['_mean']:>9.3f}\n"
           f"{'Standard Dev':<21}{stats['_std']:>9.3f}\n"
           f"{'Station Skew':<21}{stats.get('Station Skew', np.nan):>9.3f}\n"
           f"{'Regional Skew':<21}{stats.get('Regional Skew', np.nan):>9.3f}\n"
           f"{'Adopted Skew':<21}{stats['_skew']:>9.3f}\n"
           f"{'Systematic Events':<21}"
           f"{int(stats.get('Systematic Events', 0)):>9d}\n"
           f"{'Historic Period':<21}"
           f"{int(stats.get('Historic Period', 0)):>9d}\n"
           f"{'Equiv Record Length':<21}"
           f"{stats.get('Equivalent Record Length', np.nan):>9.1f}")
    ax.text(0.985, 0.03, box, transform=ax.transAxes, ha="right",
            va="bottom", family="monospace", fontsize=7,
            bbox=dict(boxstyle="square,pad=0.4", facecolor="white",
                      edgecolor="0.4", linewidth=0.6))

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Wrote {out_png}")


###############################################################################
# TABLES


def build_parameters_table(data):
    """Distribution parameters and event counts, one column per duration."""
    rows = PARAM_ROWS + DIAG_ROWS + EVENT_ROWS
    out = pd.DataFrame(index=rows, columns=list(SSP_ANALYSES), dtype=object)
    for label in SSP_ANALYSES:
        if label not in data:
            continue
        s = data[label]["stats"]
        for r in rows:
            v = s.get(r, np.nan)
            if isinstance(v, float) and np.isfinite(v):
                out.at[r, label] = (f"{int(v)}" if r in EVENT_ROWS
                                    else f"{v:.3f}")
    out.index.name = "Parameter"
    path = os.path.join(OUT_DIR, "CAS_Unreg_distribution_parameters.csv")
    out.to_csv(path)
    print(f"Wrote {path}")
    return out


def build_appendix_table(data):
    out = []
    for label in SSP_ANALYSES:
        if label not in data:
            continue
        df = data[label]["curve"]
        if df.empty:
            continue
        aep = df["pct_chance"] / 100.0
        out.append(pd.DataFrame({
            "LocationName": LOCATION_NAME,
            "Duration": label,
            "DurationDays": DURATION_DAYS[label],
            "AEP": aep,
            "VarianceLog": df["variance_log"],
            "zScore": norm.ppf(1.0 - aep),
            "Value": df["computed"],
            "UpperConf": df["conf_05"],
            "LowerConf": df["conf_95"],
            "Expected": df["expected"],
        }))
    if not out:
        return pd.DataFrame()
    table = pd.concat(out, ignore_index=True).sort_values(
        ["DurationDays", "AEP"], kind="stable")
    path = os.path.join(OUT_DIR, "CAS_Unreg_frequency_table.csv")
    table.to_csv(path, index=False)
    print(f"Wrote {len(table)} rows -> {path}")
    return table


def report_comparison(data):
    targets = np.array([0.5, 0.1, 0.02, 0.01, 0.002])
    print(f"\n{'Duration':<9}{'AEP':>8}{'2009':>11}{'2026':>11}"
          f"{'diff':>10}{'pct':>8}")
    for label in SSP_ANALYSES:
        if label not in data:
            continue
        s = data[label]["stats"]
        a = lp3_flows(*STATS_2009[label], targets)
        b = lp3_flows(s["_mean"], s["_std"], s["_skew"], targets)
        for tp, x, y in zip(targets, a, b):
            print(f"{label:<9}{tp:>8.3f}{x:>11,.0f}{y:>11,.0f}"
                  f"{y - x:>10,.0f}{100 * (y - x) / x:>7.1f}%")


###############################################################################
# MAIN


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Reading SSP results from {SSP_DIR}\n")
    data = load_all()
    if not data:
        print("No SSP reports parsed; nothing to do.")
        return
    print()

    make_summary_plot(
        data, False, os.path.join(OUT_DIR, "CAS_Unreg_2026_frequency.png"),
        "Cowlitz River at Castle Rock (14243000)\n"
        "2026 Unregulated Frequency Curves")
    make_summary_plot(
        data, True, os.path.join(OUT_DIR, "CAS_Unreg_2009_vs_2026.png"),
        "Cowlitz River at Castle Rock (14243000)\n"
        "2009 Hydrology Restudy vs 2026 Unregulated Frequency Curves")
    for label in SSP_ANALYSES:
        if label in data:
            make_duration_plot(
                data, label,
                os.path.join(OUT_DIR, f"CAS_Unreg_freq_{label}.png"))

    print()
    params = build_parameters_table(data)
    build_appendix_table(data)
    print()
    print(params.to_string())
    report_comparison(data)


if __name__ == "__main__":
    main()
