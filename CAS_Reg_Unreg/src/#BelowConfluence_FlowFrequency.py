#BelowConfluence_FlowFrequency.py
# -*- coding: utf-8 -*-
"""
FINAL regulated peak flow-frequency curves for the Cowlitz at four locations
between the Castle Rock gage and the Coweeman confluence.

THE METHOD, IN ONE LINE

    Q(AEP, site) = CAS_reg(AEP)
                   + CAS_unreg(AEP) x (DA_site - DA_gage)/DA_gage x LAG_FACTOR

    The middle term is the LOCAL contribution: the incremental drainage area
    between the gage and the site, scaled off the unregulated curve by simple
    area ratio, then reduced because the local tributaries do not crest at the
    same moment the regulated mainstem does.

WHY THE LOCAL SCALES OFF THE *UNREGULATED* CURVE
    The Coweeman, Ostrander and Arkansas all enter below Mossyrock and Riffe.
    They respond to the storm, not to the release, so the storm-scale curve --
    the UNREGULATED one -- is what governs them. Scaling them off the
    regulated curve instead would shrink them in proportion to a reservoir
    upstream of them, about 18% low at the 1,000-yr.

    Pairing CAS_unreg(AEP) with CAS_reg(AEP) is not a coincidence assumption.
    It is the same river in the same event, and the regulated curve was
    DERIVED from the unregulated one by routing in #Unreg_Reg_Curve.py.
    Nothing here assumes a tributary is simultaneously at its own AEP.

INCREMENTAL AREA, NOT TRIBUTARY AREA
    Each site's local term uses the FULL incremental area between it and the
    gage, not just the named tributary's basin. Between Castle Rock and the
    Coweeman confluence that is 247 sq mi, against 197.5 sq mi for the three
    named basins alone -- the difference is ungaged local drainage along both
    banks, which contributes whether or not it has a name.

    site                      DA      incremental   local as % of gage DA
    Castle Rock gage        2229              --           --
    below Arkansas Creek    2278              49          2.2%
    below Ostrander Creek   2335             106          4.8%
    below Coweeman River    2476             247         11.1%

WHY A PLAIN AREA RATIO
    Tested three ways and it holds up at the magnitudes that matter:
      - USGS 14245000 annual peaks, WY1950-1996, paired same-storm against the
        routed unregulated Castle Rock peak: the Coweeman runs about 1.5x its
        area share at common events but converges monotonically toward 1.0x as
        events grow. WY1996 -- a 212,245 cfs unregulated event, 92% of the
        1,000-yr flow -- sits at 1.04x. See #Coweeman_HistoricPeakRatio.py.
      - Ecology gage 26C075 storm events WY2007-WY2019: 1.11x in the >60k bin,
        and that figure is a LOWER bound because the gage's rating ceiling
        censors its largest events. See #Coweeman_Proportion.py.
      - PRISM basin precipitation over the StreamStats delineations, which
        tests the equal-depth assumption the area ratio rests on directly.
        See #Coweeman_PRISM_PrecipRatio.py.
    At the design magnitude the ratio is at area proportionality, so no
    unit-runoff uplift is applied.

WHY LAG_FACTOR = 0.80
    Measured ratio of Coweeman flow at the REGULATED Castle Rock crest to its
    own event peak, over 78 storm events:

        20-40k   n=52   median 0.806
        40-60k   n=19   median 0.781
        >60k     n= 7   median 0.413   (mean 0.511, range 0.36-0.81)
        ALL      n=78   median 0.789

    0.80 is where the two well-sampled bins group, and it is the value the
    2009 study used, so it carries precedent as well as data. It is very
    likely CONSERVATIVE: a higher lag factor adds more local flow, and the
    seven largest events in this record sit well below it.

    Not adopting the tail's 0.41 is deliberate. That figure rests on seven
    storms, its mechanism is lead time rather than magnitude (lag vs |lead|
    rho -0.55 p<0.0001; lead vs event size rho -0.04 p=0.73), and the three
    events whose Coweeman crest exceeded the gage rating are missing from it
    entirely -- and those are the events where the tributary was largest.
    See #Coweeman_LagFactor_Evidence.py for every event behind both numbers.

INPUT
    ../output/regulated_frequency_inferred.csv   (from #Unreg_Reg_Curve.py)

OUTPUT
    ../output/below_confluence_frequency.csv
    ../output/diagnostics/below_confluence_frequency.png
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
REG_CSV = r"../output/regulated_frequency_inferred.csv"
OUT_CSV = r"../output/below_confluence_frequency.csv"
PLOT_PNG = r"../output/diagnostics/below_confluence_frequency.png"
# One figure per site for the memo, keyed by the same names as LOCATIONS.
SITE_PNG = r"../output/diagnostics/freq_%s.png"

GAGE_DA = 2238.0                 # Cowlitz at the Castle Rock gage, sq mi

# Ordered downstream. Each entry is the TOTAL drainage area at that point;
# the local term is the difference from GAGE_DA, so the areas do the work and
# no tributary basin is named twice.
LOCATIONS = [
    ("Castle Rock gage",      2238.0),
    ("below Arkansas Creek",  2278.0),
    ("below Ostrander Creek", 2335.0),
    ("below Coweeman River",  2476.0),
]

LAG_FACTOR = 0.80                # flat. See WHY LAG_FACTOR = 0.80 above.

TARGET_AEP = 0.001
LAG_SENSITIVITY = (0.41, 0.60, 0.80, 1.00)

COLORS = ["#1a4f8a", "#4c8c4a", "#d99b30", "#b7410e"]

# ----------------------------------------------------------------------------


# Frequency plots use the standard-normal spacing of probability paper, but
# neither axis is LABELLED in z -- the reader gets return interval below and
# annual exceedance probability above, the same pair on the same positions.
AEP_TICKS = [0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01,
             0.005, 0.002, 0.001]
# Fixed to the Section 5 figures so the two sets of curves read the same way
# -- DQC comment on Figure 8-1.
AEP_LIMITS = (0.99, 0.001)
MIN_LABELLED_RETURN_INTERVAL = 2.0
FLOW_LIMITS = (10000.0, 400000.0)


def prob_axis(ax, label_bottom=True):
    """Return interval on the bottom, AEP percent on top.

    Return intervals shorter than two years are left unlabelled -- they are
    not conventionally shown and they crowd the left end of the axis.
    """
    zt = stats.norm.ppf(1 - np.array(AEP_TICKS))
    ri = [a for a in AEP_TICKS if 1 / a >= MIN_LABELLED_RETURN_INTERVAL]
    ax.set_xticks(stats.norm.ppf(1 - np.array(ri)))
    ax.set_xticklabels(["%g" % (1 / a) for a in ri], rotation=45, fontsize=8)
    if label_bottom:
        ax.set_xlabel("Return interval (years)")
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(zt)
    top.set_xticklabels(["%g" % (a * 100) for a in AEP_TICKS], rotation=45,
                        fontsize=8)
    top.set_xlabel("Annual exceedance probability (%)")
    return top


def local_fraction(site_da):
    """Incremental area as a fraction of the gage's own drainage area."""
    return (site_da - GAGE_DA) / GAGE_DA


def build(reg):
    out = pd.DataFrame({"AEP": reg["AEP"]})
    out["cowlitz_unreg_cfs"] = reg["unreg_computed_cfs"]
    out["cowlitz_reg_cfs"] = reg["reg_inferred_cfs"]

    for name, site_da in LOCATIONS:
        key = name.lower().replace(" ", "_")
        frac = local_fraction(site_da)
        local = reg["unreg_computed_cfs"].values * frac * LAG_FACTOR
        out["%s_local_cfs" % key] = local
        out["%s_cfs" % key] = reg["reg_inferred_cfs"].values + local
        # Band: the Castle Rock regulated band, translated by the local
        # contribution. The local term carries uncertainty of its own, but it
        # is a few percent of a quantity whose own 95% band is already wider
        # than the flow itself -- adding it would not be visible on the plot
        # and would imply a precision this method does not have. What is
        # shown is therefore the GAGE's uncertainty carried downstream, and
        # the memo says so.
        out["%s_lower_cfs" % key] = reg["reg_lower_95pct_cfs"].values + local
        out["%s_upper_cfs" % key] = reg["reg_upper_95pct_cfs"].values + local
    return out


def report(out):
    print("=" * 78)
    print("SITES")
    print("=" * 78)
    print("%-24s %8s %12s %10s" % ("site", "DA sq mi", "incremental", "local/gage"))
    for name, site_da in LOCATIONS:
        print("%-24s %8.0f %12s %9.1f%%"
              % (name, site_da,
                 "--" if site_da == GAGE_DA else "%.0f" % (site_da - GAGE_DA),
                 100 * local_fraction(site_da)))
    print("\n   lag factor %.2f applied to every local term" % LAG_FACTOR)

    print("\n" + "=" * 78)
    print("REGULATED PEAK FLOW FREQUENCY (cfs)")
    print("=" * 78)
    head = "%8s %11s" % ("AEP", "CAS unreg")
    for name, _ in LOCATIONS:
        head += " %14s" % name.replace("below ", "< ")[:14]
    print(head)
    for _, r in out.iterrows():
        line = "%8.4f %11s" % (r["AEP"], format(int(r["cowlitz_unreg_cfs"]), ","))
        for name, _ in LOCATIONS:
            key = name.lower().replace(" ", "_")
            line += " %14s" % format(int(r["%s_cfs" % key]), ",")
        print(line)

    row = out.iloc[(out["AEP"] - TARGET_AEP).abs().idxmin()]
    print("\n" + "=" * 78)
    print("AT AEP=%.4f  (1-in-%d)" % (row["AEP"], round(1 / row["AEP"])))
    print("=" * 78)
    print("   Castle Rock unregulated %11s cfs   (drives every local term)"
          % format(int(row["cowlitz_unreg_cfs"]), ","))
    print("%-24s %12s %11s %11s   %s"
          % ("site", "local cfs", "TOTAL", "over gage", "95% band"))
    base = row["castle_rock_gage_cfs"]
    for name, _ in LOCATIONS:
        key = name.lower().replace(" ", "_")
        total = row["%s_cfs" % key]
        print("%-24s %12s %11s %10.1f%%   %s to %s"
              % (name, format(int(row["%s_local_cfs" % key]), ","),
                 format(int(total), ","), 100 * (total - base) / base,
                 format(int(row["%s_lower_cfs" % key]), ","),
                 format(int(row["%s_upper_cfs" % key]), ",")))

    print("\n" + "=" * 78)
    print("SENSITIVITY TO LAG_FACTOR, at the Coweeman confluence, AEP=%.4f"
          % TARGET_AEP)
    print("=" * 78)
    frac = local_fraction(LOCATIONS[-1][1])
    band = row["below_coweeman_river_upper_cfs"] - row["below_coweeman_river_lower_cfs"]
    adopted = row["below_coweeman_river_cfs"]
    for trial in LAG_SENSITIVITY:
        total = row["cowlitz_reg_cfs"] + row["cowlitz_unreg_cfs"] * frac * trial
        mark = "   <- adopted" if abs(trial - LAG_FACTOR) < 1e-9 else ""
        print("   lag %.2f -> local %8s   total %11s  (%+.1f%% vs adopted)%s"
              % (trial, format(int(row["cowlitz_unreg_cfs"] * frac * trial), ","),
                 format(int(total), ","), 100 * (total - adopted) / adopted, mark))
    print("\n   the whole 0.41-1.00 range spans %.1f%% of the adopted total, "
          "against\n   a 95%% band %s cfs wide (%.0f%% of it)."
          % (100 * (row["cowlitz_unreg_cfs"] * frac * (1.00 - 0.41)) / adopted,
             format(int(band), ","),
             100 * (row["cowlitz_unreg_cfs"] * frac * (1.00 - 0.41)) / band))


def plot(out):
    """Figure 8-1: the four regulated curves, nothing else.

    The uncertainty band and the percentage-increase panel were both dropped
    at DQC review -- the band is Castle Rock's, carried forward unchanged and
    already shown in Section 5, and the percentages are in Table 8-3.
    """
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, ax = plt.subplots(figsize=(10, 7.5))

    ax.set_yscale("log")   # BEFORE any annotation -- a pre-log get_ylim()
                           # grabs a linear autoscale limit and crushes the
                           # data into a sliver once the scale changes.
    for (name, site_da), color in zip(LOCATIONS, COLORS):
        key = name.lower().replace(" ", "_")
        ax.plot(z, out["%s_cfs" % key], color=color, lw=2.2,
                label="%s  (%.0f sq mi)" % (name, site_da))

    ax.set_ylim(FLOW_LIMITS)
    ax.set_xlim(stats.norm.ppf(1 - AEP_LIMITS[0]),
                stats.norm.ppf(1 - AEP_LIMITS[1]))
    ax.set_ylabel("Regulated peak flow (cfs)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: format(int(v), ",")))
    ax.set_title("Regulated peak flow frequency, Castle Rock gage to the "
                 "Coweeman confluence", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    prob_axis(ax, label_bottom=True)

    fig.tight_layout()
    os.makedirs(os.path.dirname(PLOT_PNG), exist_ok=True)
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


def site_plot(out, name, site_da, color):
    """One site, one figure. Castle Rock also carries the unregulated curve
    because that is the pair the reader needs there; downstream the
    unregulated curve is an input to the local term, not a result, so it
    would only invite the wrong comparison."""
    key = name.lower().replace(" ", "_")
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, ax = plt.subplots(figsize=(9, 6.4))
    ax.set_yscale("log")        # before any annotation -- see plot()
    ax.fill_between(z, out["%s_lower_cfs" % key], out["%s_upper_cfs" % key],
                    color=color, alpha=0.14,
                    label="95% confidence band")
    if site_da == GAGE_DA:
        ax.plot(z, out["cowlitz_unreg_cfs"], color="#7aa9d0", lw=2, ls="--",
                label="Unregulated")
    ax.plot(z, out["%s_cfs" % key], color=color, lw=2.6, label="Regulated")
    ax.axvline(stats.norm.ppf(1 - TARGET_AEP), color="gray", lw=1, ls=":")
    ax.set_ylabel("Peak flow (cfs)")
    sub = ("Cowlitz at Castle Rock: Regulated and unregulated peak flow frequency"
           if site_da == GAGE_DA else
           "Cowlitz %s regulated peak flow frequency" % name)
    ax.set_title(sub)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    prob_axis(ax)

    fig.tight_layout()
    path = SITE_PNG % key
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_site_tables(out):
    """One CSV per site: the table that goes in the memo, nothing else."""
    paths = []
    for name, _ in LOCATIONS:
        key = name.lower().replace(" ", "_")
        frame = pd.DataFrame({
            "AEP": out["AEP"],
            "regulated_cfs": out["%s_cfs" % key].round(0),
            "lower_95pct_cfs": out["%s_lower_cfs" % key].round(0),
            "upper_95pct_cfs": out["%s_upper_cfs" % key].round(0),
        })
        if key == "castle_rock_gage":
            frame.insert(1, "unregulated_cfs",
                         out["cowlitz_unreg_cfs"].round(0))
        else:
            frame.insert(1, "local_cfs", out["%s_local_cfs" % key].round(0))
        path = r"../output/freq_table_%s.csv" % key
        frame.to_csv(path, index=False)
        paths.append(path)
    return paths


def main():
    reg = pd.read_csv(REG_CSV).sort_values("AEP", ascending=False).reset_index(drop=True)
    out = build(reg)
    out.to_csv(OUT_CSV, index=False)
    report(out)
    plot(out)
    for (name, site_da), color in zip(LOCATIONS, COLORS):
        print("Wrote", site_plot(out, name, site_da, color))
    for path in write_site_tables(out):
        print("Wrote", path)
    print("Wrote", OUT_CSV)


if __name__ == "__main__":
    main()
