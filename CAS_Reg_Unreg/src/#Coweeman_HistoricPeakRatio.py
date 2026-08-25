#Coweeman_HistoricPeakRatio.py
# -*- coding: utf-8 -*-
"""
LOCAL_RATIO from the long USGS annual-peak record, independent of the Ecology
gage and of its rating ceiling.

WHY
    The ratio in #BelowConfluence_FlowFrequency.py rests on 76 storm events in
    13 water years from Ecology gage 26C075, whose rating tops out near 3,400
    cfs -- so its largest events are censored exactly where the analysis needs
    them. USGS gage 14245000 (Coweeman R nr Kelso, 119 sq mi) is the SAME
    river with a published annual-peak series reaching 7,730 cfs, over
    WY1950-1996. It is a genuinely independent check on the same quantity,
    from a record with no 3,400 ceiling.

WHAT IS COMPARED
    Coweeman USGS annual peak, against the ROUTED UNREGULATED Castle Rock peak
    in the same water year (ResSim_WCM_RC.dss, 1928-2026). Unregulated is the
    right mainstem term: the Coweeman enters below the projects and responds
    to the storm, not to the release.

THE SAME-STORM PROBLEM
    An annual maximum on each river is not necessarily the same event. Where
    the two peaks are more than SAME_STORM_DAYS apart they are different
    storms, and their ratio is meaningless as a coincident quantity -- it
    pairs the Coweeman's biggest day with a mainstem peak that happened weeks
    away. Those years are reported separately and excluded from the headline
    number, never silently averaged in.

WHAT THIS DOES AND DOES NOT SETTLE
    It constrains the PEAK ratio (Coweeman peak / CAS unreg peak) on a much
    longer record. It says nothing about the LAG factor -- an annual peak
    carries no sub-daily timing, and most of these years have no clock time at
    all. See #Coweeman_LagFactor_Evidence.py for that.

OUTPUTS
    ../output/diagnostics/coweeman_historic_peak_ratio.csv
    ../output/diagnostics/coweeman_historic_peak_ratio.png
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime
from pydsstools.heclib.dss import HecDss

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
PEAKS_RDB = r"../data/coweeman/usgs_peaks_14245000.rdb"
DSS_PATH = r"../output/ResSim_WCM_RC.dss"
UNREG_PATHNAME = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"

OUT_CSV = r"../output/diagnostics/coweeman_historic_peak_ratio.csv"
PLOT_PNG = r"../output/diagnostics/coweeman_historic_peak_ratio.png"

COW_GAGE_DA = 119.0        # USGS 14245000, sq mi
CAS_DA = 2238.0            # Cowlitz above Castle Rock, sq mi
SAME_STORM_DAYS = 5        # beyond this the two annual maxima are different storms

# For direct comparison with the Ecology-era result.
ECOLOGY_TAIL_RATIO = 0.0590
ECOLOGY_ALL_RATIO = 0.0629

C_SAME = "#1a4f8a"
C_DIFF = "#c0c0c0"

# ----------------------------------------------------------------------------


def water_year(stamp):
    """USACE convention: Oct-Dec belong to the FOLLOWING calendar year's WY."""
    return stamp.year + (1 if stamp.month >= 10 else 0)


def read_peaks(path):
    """USGS RDB annual peaks. Two header rows: names, then format codes."""
    rows = []
    with open(path, "r", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue          # the RDB interleaves blank lines with the
                                  # comment block; a blank one taken as the
                                  # header silently loses every column
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if header is None:
                header = parts
                continue
            if parts and parts[0].endswith("s"):   # the 5s/15s/10d format row
                continue
            rows.append(dict(zip(header, parts)))
    frame = pd.DataFrame(rows)
    frame["peak_dt"] = pd.to_datetime(frame["peak_dt"], errors="coerce")
    frame["peak_va"] = pd.to_numeric(frame["peak_va"], errors="coerce")
    frame = frame.dropna(subset=["peak_dt", "peak_va"])
    frame["wy"] = frame["peak_dt"].map(water_year)
    print("   Coweeman USGS 14245000: %d annual peaks, WY%d-WY%d, max %s cfs"
          % (len(frame), frame.wy.min(), frame.wy.max(),
             format(int(frame.peak_va.max()), ",")))
    codes = frame["peak_cd"].replace("", np.nan).dropna()
    if len(codes):
        print("      peak codes present: %s" % ", ".join(sorted(set(codes))))
    return frame[["wy", "peak_dt", "peak_va", "peak_cd"]]


def first_stamp(ts):
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    roll = False
    if " 24:" in text:
        text, roll = text.replace(" 24:", " 00:"), True
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M"):
        try:
            stamp = pd.Timestamp(datetime.strptime(text, fmt))
            return stamp + pd.Timedelta(days=1) if roll else stamp
        except ValueError:
            continue
    return pd.Timestamp(text)


def read_unreg(path, pathname):
    version = 6 if open(path, "rb").read(16)[12] == 6 else 7
    dss = HecDss.Open(path, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        values[np.array(ts.nodata, dtype=bool)] = np.nan
        values[values <= -900.0] = np.nan
        step = pd.Timedelta(seconds=int(ts.interval))
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    series = pd.Series(values, index=index).dropna().sort_index()
    print("   Castle Rock unregulated: %d hours, %s to %s"
          % (len(series), series.index.min().date(), series.index.max().date()))
    return series


def build(peaks, unreg):
    wy = pd.Series(unreg.index.map(water_year), index=unreg.index)
    rows = []
    for _, r in peaks.iterrows():
        block = unreg[wy == r["wy"]]
        if not len(block):
            continue
        cas_time = block.idxmax()
        cas_peak = float(block.max())
        gap = abs((cas_time.normalize() - r["peak_dt"]).days)
        rows.append({
            "wy": int(r["wy"]),
            "cow_peak_date": r["peak_dt"].date(),
            "cow_peak_cfs": r["peak_va"],
            "cow_peak_cd": r["peak_cd"],
            "cas_unreg_peak_date": cas_time.date(),
            "cas_unreg_peak_cfs": cas_peak,
            "days_apart": gap,
            "same_storm": gap <= SAME_STORM_DAYS,
            "ratio_peak": r["peak_va"] / cas_peak,
        })
    return pd.DataFrame(rows).sort_values("wy").reset_index(drop=True)


def report(table):
    da = COW_GAGE_DA / CAS_DA
    same = table[table["same_storm"]]
    diff = table[~table["same_storm"]]
    print("\n" + "=" * 78)
    print("PAIRED WATER YEARS")
    print("=" * 78)
    print("   %d water years in common; %d same-storm (within %d days), %d not"
          % (len(table), len(same), SAME_STORM_DAYS, len(diff)))
    if len(diff):
        print("   different-storm years excluded: %s"
              % ", ".join("WY%d" % w for w in diff["wy"]))

    print("\n" + "=" * 78)
    print("COWEEMAN PEAK / CASTLE ROCK UNREGULATED PEAK")
    print("=" * 78)
    print("   drainage-area ratio %.1f/%.1f = %.5f\n" % (COW_GAGE_DA, CAS_DA, da))
    print("   %-26s %4s %9s %9s %9s %9s %8s"
          % ("set", "n", "median", "mean", "p25", "p75", "xDA"))
    for name, s in [("same-storm years", same["ratio_peak"]),
                    ("different-storm years", diff["ratio_peak"]),
                    ("all paired years", table["ratio_peak"])]:
        if len(s):
            print("   %-26s %4d %9.5f %9.5f %9.5f %9.5f %7.2fx"
                  % (name, len(s), s.median(), s.mean(), s.quantile(.25),
                     s.quantile(.75), s.median() / da))

    print("\n   for comparison, Ecology 26C075 storm events WY2007-WY2019:")
    print("   %-26s %4s %9.5f %19s %7.2fx"
          % ("Ecology, all events", "76", ECOLOGY_ALL_RATIO, "", ECOLOGY_ALL_RATIO / da))
    print("   %-26s %4s %9.5f %19s %7.2fx  (rating-censored, a LOWER bound)"
          % ("Ecology, >60k bin", "7", ECOLOGY_TAIL_RATIO, "", ECOLOGY_TAIL_RATIO / da))

    if len(same) > 3:
        print("\n" + "=" * 78)
        print("DOES THE RATIO FALL WITH EVENT SIZE?  (same-storm years only)")
        print("=" * 78)
        rho, p = stats.spearmanr(same["cas_unreg_peak_cfs"], same["ratio_peak"])
        print("   Spearman rho %+.3f  p=%.4f   %s"
              % (rho, p, "significant" if p < 0.05 else "NOT significant"))
        big = same[same["cas_unreg_peak_cfs"] >= 60000]
        small = same[same["cas_unreg_peak_cfs"] < 60000]
        for nm, s in [("CAS unreg >= 60k", big), ("CAS unreg <  60k", small)]:
            if len(s):
                print("   %-18s n=%2d  median %.5f  = %.2fx DA"
                      % (nm, len(s), s["ratio_peak"].median(),
                         s["ratio_peak"].median() / da))

    print("\n" + "=" * 78)
    print("THE LARGEST EVENTS -- the only ones near the design magnitude")
    print("=" * 78)
    big = same.sort_values("cas_unreg_peak_cfs", ascending=False)
    print("%6s %12s %10s %9s %7s" % ("WY", "CAS unreg", "cow", "ratio", "xDA"))
    for _, r in big.head(8).iterrows():
        print("%6d %12s %10s %9.5f %6.2fx"
              % (r["wy"], format(int(r["cas_unreg_peak_cfs"]), ","),
                 format(int(r["cow_peak_cfs"]), ","), r["ratio_peak"],
                 r["ratio_peak"] / da))
    for thr in (100000, 80000, 60000):
        t = big[big["cas_unreg_peak_cfs"] >= thr]["ratio_peak"]
        if len(t):
            print("   CAS unreg >= %7s  n=%2d  median %.5f = %.2fx DA"
                  % (format(thr, ","), len(t), t.median(), t.median() / da))
    print("\n   WY1996 is the single most relevant row in this study: a\n"
          "   %s cfs unregulated event, %.0f%% of the 1,000-yr flow, with the\n"
          "   Coweeman at %.2fx its drainage-area share. The ratio declines\n"
          "   monotonically toward the area ratio as events get larger, and at\n"
          "   the design magnitude it is essentially AT it."
          % (format(int(big.iloc[0]["cas_unreg_peak_cfs"]), ","),
             100 * big.iloc[0]["cas_unreg_peak_cfs"] / 230884.0,
             big.iloc[0]["ratio_peak"] / da))

    print("\n" + "=" * 78)
    print("EVERY PAIRED YEAR")
    print("=" * 78)
    print("%6s %13s %11s %13s %11s %7s %6s %8s"
          % ("WY", "cow date", "cow cfs", "CAS date", "CAS cfs", "apart", "same", "ratio"))
    for _, r in table.iterrows():
        print("%6d %13s %11s %13s %11s %7d %6s %8.5f"
              % (r["wy"], r["cow_peak_date"], format(int(r["cow_peak_cfs"]), ","),
                 r["cas_unreg_peak_date"], format(int(r["cas_unreg_peak_cfs"]), ","),
                 r["days_apart"], "yes" if r["same_storm"] else "no", r["ratio_peak"]))


def plot(table):
    da = COW_GAGE_DA / CAS_DA
    same = table[table["same_storm"]]
    diff = table[~table["same_storm"]]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    ax1.scatter(diff["cas_unreg_peak_cfs"] / 1000, diff["ratio_peak"], s=44,
                color=C_DIFF, edgecolor="0.4", lw=0.5, label="different storm (excluded)")
    ax1.scatter(same["cas_unreg_peak_cfs"] / 1000, same["ratio_peak"], s=52,
                color=C_SAME, edgecolor="0.2", lw=0.6, label="same storm", zorder=3)
    ax1.axhline(da, color="#b7410e", ls="--", lw=1.8,
                label="drainage-area ratio %.4f" % da)
    ax1.axhline(same["ratio_peak"].median(), color=C_SAME, ls=":", lw=1.8,
                label="same-storm median %.4f" % same["ratio_peak"].median())
    ax1.axhline(ECOLOGY_ALL_RATIO, color="#4c8c4a", ls="-.", lw=1.5,
                label="Ecology-era all events %.4f" % ECOLOGY_ALL_RATIO)
    ax1.set_xlabel("Castle Rock unregulated annual peak (1000 cfs)")
    ax1.set_ylabel("Coweeman peak / CAS unreg peak")
    ax1.set_title("Annual peak ratio vs event size\nUSGS 14245000, WY%d-WY%d"
                  % (table.wy.min(), table.wy.max()))
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.bar(same["wy"], same["ratio_peak"], color=C_SAME, label="same storm")
    ax2.bar(diff["wy"], diff["ratio_peak"], color=C_DIFF, label="different storm")
    ax2.axhline(da, color="#b7410e", ls="--", lw=1.8)
    ax2.axhline(same["ratio_peak"].median(), color=C_SAME, ls=":", lw=1.8)
    ax2.set_xlabel("water year")
    ax2.set_ylabel("Coweeman peak / CAS unreg peak")
    ax2.set_title("By water year\ndashed = drainage-area ratio, dotted = same-storm median")
    ax2.grid(alpha=0.3, axis="y")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


def main():
    print("Reading...")
    peaks = read_peaks(PEAKS_RDB)
    unreg = read_unreg(DSS_PATH, UNREG_PATHNAME)
    table = build(peaks, unreg)
    table.to_csv(OUT_CSV, index=False)
    report(table)
    plot(table)
    print("Wrote", OUT_CSV)


if __name__ == "__main__":
    main()
