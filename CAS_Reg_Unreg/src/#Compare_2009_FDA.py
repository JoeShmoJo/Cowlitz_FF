#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#Compare_2009_FDA.py

The 2026 regulated frequency curves against the 2009 Restudy, with both
studies' uncertainty, at every location the 2009 study exported from HEC-FDA.

    output/compare_2009_<location>.csv              appendix tables
    output/diagnostics/compare_2009_<location>.png  appendix figures

The Castle Rock figure is the one that goes in the main report. The other two
are appendix only. All three are drawn the same way so they can be read
against each other.

WHY THE 2009 CURVE COMES FROM HEC-FDA AND NOT FROM TABLE B-5
------------------------------------------------------------
Table B-5 of the 2009 Restudy publishes a regulated curve for Castle Rock
without uncertainty. The HEC-FDA 2.0.2 Cowlitz study carries a different curve
for the same site, on the same AEP grid, and it is the one with uncertainty
attached. The author adopted the FDA curve. Table B-5 is preserved unused at
data/fda_2009/curve_2009_tableB5_superseded.csv.

The two differ by -49% at the 1.01 year and +23% near the 50 year and agree
only at the two rarest ordinates. Ranked against the regulated records, mean
absolute error in log10:

                     vs OBSERVED (n=51)    vs ADJUSTED, FRM only (n=41)
  2026 regulated        0.0748  18.8%          0.0593  14.6%
  2009 Table B-5        0.0819  20.8%          0.0462  11.2%
  2009 HEC-FDA          0.0267   6.3%          0.1040  27.1%

Table B-5 tracks the adjusted, flood risk management only record. The FDA
curve tracks the observed record, which still carries the incidental reservoir
drafting that the adjustment in Section 5.2 exists to remove. Worth carrying
forward: the comparison in these figures is therefore between a 2026 curve on
an FRM only basis and a 2009 curve that is not, and the gap at the frequent end
is partly that difference in basis rather than a difference between studies.

THE BANDS ARE PUT ON A COMMON FOOTING
    FDA reports its bounds at 2.5 and 97.5 percent, a 95 percent interval. The
    memo is 90 percent throughout. Drawing both as reported would put two
    interval widths on one figure, which is the confusion the DQC review asked
    to have removed, so the FDA bounds are rescaled. They are log-symmetric
    about their curve to within 0.003 dex, so each side is reduced to its own
    sigma and re-expressed at 90 percent with nothing lost.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
FDA_DIR = r"../data/fda_2009"
OUT_DIR = r"../output"
DIAG_DIR = r"../output/diagnostics"

# location key -> (pretty name, FDA export, 2026 source, 2026 flow column)
# Castle Rock's 2026 curve comes from the transform output; the downstream
# locations come from the confluence model, which adds the local contribution.
LOCATIONS = [
    ("castle_rock_gage", "Castle Rock gage",
     "fda_2009_castle_rock_gage.csv",
     "regulated_frequency_inferred.csv", "reg_inferred_cfs",
     "reg_lower_90pct_cfs", "reg_upper_90pct_cfs"),
    ("below_arkansas_creek", "Below Arkansas Creek",
     "fda_2009_below_arkansas_creek.csv",
     "freq_table_below_arkansas_creek.csv", "regulated_cfs",
     "lower_90pct_cfs", "upper_90pct_cfs"),
    ("below_ostrander_creek", "Below Ostrander Creek",
     "fda_2009_below_ostrander_creek.csv",
     "freq_table_below_ostrander_creek.csv", "regulated_cfs",
     "lower_90pct_cfs", "upper_90pct_cfs"),
]

# The memo's interval. Matches UNCERTAINTY_CONF_LEVEL in #Unreg_Reg_Curve.py
# and must stay in step with it.
CONF_LEVEL = 0.90
# What FDA reported, so the rescale knows what it is undoing.
FDA_CONF_LEVEL = 0.95

# Only these AEPs go in the appendix tables. The full grid is in the CSVs.
TABLE_AEPS = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005,
              0.0001]

AEP_TICKS = [0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005,
             0.002, 0.001]
AEP_LIMITS = (0.99, 0.001)
MIN_LABELLED_RETURN_INTERVAL = 2.0
FLOW_LIMITS = (10000.0, 400000.0)
FIG_SIZE = (10.5, 8.2)

C_2026 = "#c0392b"
C_2009 = "#1a7a5e"

# ----------------------------------------------------------------------------


def probability_axis(ax):
    """SSP-style normal-probability axis, AEP decreasing to the right.

    Return interval along the bottom and annual exceedance probability across
    the top, the same arrangement #Unreg_Reg_Curve.py uses, so these figures
    and Figure 6-2 can be read against each other. Return intervals shorter
    than two years are left unlabelled, per the DQC comment on Figure 6-1.
    """
    ax.set_xlim(stats.norm.ppf(1.0 - AEP_LIMITS[0]),
                stats.norm.ppf(1.0 - AEP_LIMITS[1]))
    ri = [t for t in AEP_TICKS if 1.0 / t >= MIN_LABELLED_RETURN_INTERVAL]
    ax.xaxis.set_major_locator(
        FixedLocator(stats.norm.ppf(1.0 - np.array(ri))))
    ax.set_xticklabels(["%g" % (1.0 / t) for t in ri], rotation=45, fontsize=8)
    ax.set_xlabel("Return interval (years)")

    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.xaxis.set_major_locator(
        FixedLocator(stats.norm.ppf(1.0 - np.array(AEP_TICKS))))
    top.set_xticklabels(["%g" % (t * 100) for t in AEP_TICKS], rotation=45,
                        fontsize=8)
    top.set_xlabel("Annual exceedance probability (%)")


def read_fda(name):
    """One FDA export, with its bounds rescaled to the memo's interval.

    Each side keeps its own sigma, so the asymmetry FDA reported survives the
    rescale instead of being averaged away.
    """
    path = os.path.join(FDA_DIR, name)
    if not os.path.exists(path):
        raise SystemExit("missing FDA export: %s" % path)
    f = pd.read_csv(path, comment="#")
    z_fda = stats.norm.ppf(0.5 + FDA_CONF_LEVEL / 2.0)
    z_ours = stats.norm.ppf(0.5 + CONF_LEVEL / 2.0)
    sig_lo = np.log10(f["cfs"] / f["lo_2p5"]) / z_fda
    sig_hi = np.log10(f["hi_97p5"] / f["cfs"]) / z_fda
    asym = np.abs(sig_hi - sig_lo).max()
    out = pd.DataFrame({
        "AEP": f["aep"].values,
        "cfs_2009": f["cfs"].values,
        "lo_2009": (f["cfs"] * 10.0 ** (-z_ours * sig_lo)).values,
        "hi_2009": (f["cfs"] * 10.0 ** (z_ours * sig_hi)).values,
    })
    out.attrs["sigma_range"] = (sig_lo.min(), sig_hi.max())
    out.attrs["asymmetry_dex"] = asym
    return out


def read_2026(csv_name, flow_col, lo_col, hi_col):
    path = os.path.join(OUT_DIR, csv_name)
    if not os.path.exists(path):
        raise SystemExit("missing 2026 curve: %s" % path)
    d = pd.read_csv(path)
    return pd.DataFrame({"AEP": d["AEP"].values,
                         "cfs_2026": d[flow_col].values,
                         "lo_2026": d[lo_col].values,
                         "hi_2026": d[hi_col].values})


def at_aep(aep_known, flow_known, aep_want):
    """Interpolate in z against log flow, the space the curves are drawn in.

    Outside the tabulated range the result is NaN rather than an extrapolation,
    because neither study's curve should be read past its own last ordinate.
    """
    z = stats.norm.ppf(1.0 - np.asarray(aep_known, dtype=float))
    q = np.log10(np.asarray(flow_known, dtype=float))
    order = np.argsort(z)
    zw = stats.norm.ppf(1.0 - np.asarray(aep_want, dtype=float))
    return 10.0 ** np.interp(zw, z[order], q[order], left=np.nan,
                             right=np.nan)


def build_table(fda, cur):
    """Both studies on the 2026 AEP grid, joined for the appendix."""
    out = cur.copy()
    for col in ("cfs_2009", "lo_2009", "hi_2009"):
        out[col] = at_aep(fda["AEP"], fda[col], out["AEP"])
    out["diff_cfs"] = out["cfs_2026"] - out["cfs_2009"]
    out["diff_pct"] = 100.0 * (out["cfs_2026"] / out["cfs_2009"] - 1.0)
    # Do the two bands overlap at all? Where they do not, the studies disagree
    # by more than either one's own uncertainty admits, which is the thing a
    # reader of these figures most wants to know.
    out["bands_overlap"] = ((out["lo_2026"] <= out["hi_2009"]) &
                            (out["lo_2009"] <= out["hi_2026"]))
    return out


def plot_location(key, pretty, table, stem):
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    pct = int(round(100 * CONF_LEVEL))

    z26 = stats.norm.ppf(1.0 - table["AEP"].values)
    ax.fill_between(z26, table["lo_2026"], table["hi_2026"], color=C_2026,
                    alpha=0.17, lw=0, zorder=1,
                    label="2026, %d%%" % pct)
    d9 = table.dropna(subset=["cfs_2009"])
    z09 = stats.norm.ppf(1.0 - d9["AEP"].values)
    ax.fill_between(z09, d9["lo_2009"], d9["hi_2009"], color=C_2009,
                    alpha=0.13, lw=0, zorder=2,
                    label="2009 study, %d%%" % pct)
    ax.plot(z09, d9["lo_2009"], color=C_2009, lw=0.9, ls=":", zorder=3)
    ax.plot(z09, d9["hi_2009"], color=C_2009, lw=0.9, ls=":", zorder=3)

    ax.plot(z26, table["cfs_2026"], color=C_2026, lw=2.6, zorder=5,
            label="2026 regulated")
    ax.plot(z09, d9["cfs_2009"], color=C_2009, lw=1.7, ls="--", zorder=4,
            label="2009 study, regulated")

    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    probability_axis(ax)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("%s, regulated peak flow frequency\n2026 study against the "
                 "2009 Restudy, both at %d percent" % (pretty, pct),
                 fontsize=12)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    fig.tight_layout()
    fig.savefig(stem, dpi=150)
    plt.close(fig)
    print("   figure", stem)


def report(pretty, table, fda):
    sub = table.dropna(subset=["cfs_2009"])
    lo, hi = fda.attrs["sigma_range"]
    print("\n%s" % pretty)
    print("   FDA sigma %.4f to %.4f dex, bound asymmetry up to %.4f dex"
          % (lo, hi, fda.attrs["asymmetry_dex"]))
    print("   %-8s %11s %11s %9s %8s  %s"
          % ("AEP", "2026", "2009", "diff", "diff %", "bands"))
    for aep in TABLE_AEPS:
        m = sub[np.isclose(sub["AEP"], aep)]
        if not len(m):
            continue
        r = m.iloc[0]
        print("   %-8g %11s %11s %9s %7.1f%%  %s"
              % (aep, format(int(round(r.cfs_2026)), ","),
                 format(int(round(r.cfs_2009)), ","),
                 format(int(round(r.diff_cfs)), ","), r.diff_pct,
                 "overlap" if r.bands_overlap else "DISJOINT"))
    print("   median difference %+.1f%%, range %+.1f%% to %+.1f%%"
          % (sub.diff_pct.median(), sub.diff_pct.min(), sub.diff_pct.max()))
    n = int((~sub.bands_overlap).sum())
    print("   the two %d%% bands are disjoint at %d of %d shared ordinates"
          % (int(round(100 * CONF_LEVEL)), n, len(sub)))


def main():
    for d in (OUT_DIR, DIAG_DIR):
        if not os.path.isdir(d):
            os.makedirs(d)
    print("2009 HEC-FDA against the 2026 curves, %d percent bands on both"
          % int(round(100 * CONF_LEVEL)))
    for key, pretty, fda_csv, cur_csv, fcol, lcol, hcol in LOCATIONS:
        fda = read_fda(fda_csv)
        cur = read_2026(cur_csv, fcol, lcol, hcol)
        table = build_table(fda, cur)
        out_csv = os.path.join(OUT_DIR, "compare_2009_%s.csv" % key)
        table.to_csv(out_csv, index=False, float_format="%.1f")
        report(pretty, table, fda)
        print("   table ", out_csv)
        plot_location(key, pretty, table,
                      os.path.join(DIAG_DIR, "compare_2009_%s.png" % key))


main()
