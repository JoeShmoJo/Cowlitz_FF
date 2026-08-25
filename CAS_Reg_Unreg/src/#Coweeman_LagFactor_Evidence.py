#Coweeman_LagFactor_Evidence.py
# -*- coding: utf-8 -*-
"""
Every piece of evidence behind the LAG FACTOR, laid out so a number can be
chosen rather than asserted.

WHAT THE LAG FACTOR IS
    The Coweeman does not crest when the regulated Cowlitz does. It peaks
    early and is already receding when the mainstem arrives, so only part of
    its own peak is present at the moment that matters:

        lag_factor = (Coweeman flow at the REGULATED Castle Rock crest)
                     / (Coweeman peak in the same event)

    #BelowConfluence_FlowFrequency.py multiplies the local peak by this. It is
    the least-evidenced number in that script, which is why this one exists.

WHAT THIS SCRIPT DOES NOT DO
    It does not pick the value. It prints the candidates, the support behind
    each, and what each costs at the 1,000-yr, and stops.

THE THREE THINGS THAT MAKE THE TAIL VALUE HARD
    1. n=7 in the >60k bin, spanning 0.36 to 0.81. Median 0.420, mean 0.511 --
       at that sample size the two statistics disagree and neither is stable.
    2. The decline with magnitude is real (Spearman rho about -0.32,
       p about 0.004 over all events) but its MECHANISM is lead time, not
       size: lag_factor vs |lead| is the stronger relationship, and lead
       itself is NOT significantly related to event size (p about 0.53). So
       the decline may be a property of which storms happened to be large in
       a 13-year sample rather than of largeness itself.
    3. Three events are missing from the tail entirely because the Coweeman
       was ABOVE ITS RATING at the crest and Ecology reported nothing. Those
       are events where the tributary was LARGEST, so every statistic here is
       biased LOW by their absence. They are listed, not hidden.

OUTPUTS
    ../output/diagnostics/coweeman_lagfactor_events.csv
    ../output/diagnostics/coweeman_lagfactor_evidence.png
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime
from pydsstools.heclib.dss import HecDss

REPO_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "Modules"))
from ecology_io import (read_ecology_cache, MISSING_CODES, TRUSTED_CODES,
                        CODE_MEANING)

# ----------------------------------------------------------------------------
# USER SETTINGS  -- event selection HELD IDENTICAL to #Coweeman_Proportion.py
# ----------------------------------------------------------------------------
CACHE_DIR = r"../data/coweeman"
DSS_PATH = r"../output/ResSim_WCM_RC.dss"
UNREG_PATHNAME = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"
REG_PATHNAME = "//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/"

OUT_CSV = r"../output/diagnostics/coweeman_lagfactor_events.csv"
PLOT_PNG = r"../output/diagnostics/coweeman_lagfactor_evidence.png"

N_EVENTS = 80
EVENT_MIN_SEPARATION_DAYS = 7
MIN_EVENT_CFS = 20000.0
EVENT_WINDOW_HOURS = 48
REG_SEARCH_WINDOW_HOURS = 72

BINS = [20000, 40000, 60000, 1e9]
LABELS = ["20-40k", "40-60k", ">60k"]

# For costing each candidate at the design event.
TARGET_AEP = 0.001
CAS_UNREG_AT_TARGET = 230884.0
CAS_REG_AT_TARGET = 189434.0
LOCAL_RATIO = 0.09796          # from #BelowConfluence_FlowFrequency.py

C_BIN = {"20-40k": "#9bb8d4", "40-60k": "#4c7fb0", ">60k": "#1a4f8a"}

# ----------------------------------------------------------------------------


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


def read_dss_series(path, pathname, label):
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
    print("   %-24s %d hours, %s to %s" % (label, len(series),
          series.index.min().date(), series.index.max().date()))
    return series


def find_events(series, n, min_separation_days, floor):
    working, events = series.copy(), []
    span = pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span) & (working.index <= when + span)] = np.nan
    return sorted(events)


def at(series, when, tol_min=30):
    """Exact-ish lookup. NOT interpolated -- interpolation across a censored
    gap invents a number the gage never reported."""
    if not len(series):
        return np.nan
    i = series.index.get_indexer([when], method="nearest")[0]
    if i < 0 or abs(series.index[i] - when) > pd.Timedelta(minutes=tol_min):
        return np.nan
    return float(series.iloc[i])


def code_at(qual, when, tol_min=30):
    if not len(qual):
        return 0
    i = qual.index.get_indexer([when], method="nearest")[0]
    if i < 0 or abs(qual.index[i] - when) > pd.Timedelta(minutes=tol_min):
        return 0
    return int(qual.iloc[i])


def build(cow, qual, unreg, reg):
    rows = []
    for when, unreg_peak in find_events(unreg, N_EVENTS,
                                        EVENT_MIN_SEPARATION_DAYS, MIN_EVENT_CFS):
        win = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
        cow_win = cow.loc[when - win:when + win]
        reg_win = reg.loc[when - pd.Timedelta(hours=REG_SEARCH_WINDOW_HOURS):
                          when + pd.Timedelta(hours=REG_SEARCH_WINDOW_HOURS)]
        if not len(reg_win) or not cow_win.notna().any():
            continue
        cow_peak_time = cow_win.idxmax()
        cow_peak = float(cow_win.max())
        reg_peak_time = reg_win.idxmax()
        cow_at_reg = at(cow, reg_peak_time)

        # Was the gage censored anywhere across the crest window? If so the
        # true peak is above what is recorded and this row is a lower bound.
        span = qual.loc[min(cow_peak_time, reg_peak_time) - pd.Timedelta(hours=6):
                        max(cow_peak_time, reg_peak_time) + pd.Timedelta(hours=6)]
        censored = bool(span.isin(MISSING_CODES).any())

        rows.append({
            "event_time": when,
            "cas_unreg_peak_cfs": unreg_peak,
            "cas_reg_peak_cfs": float(reg_win.max()),
            "reg_minus_unreg_hours": (reg_peak_time - when).total_seconds() / 3600.0,
            "cow_peak_cfs": cow_peak,
            "cow_peak_time": cow_peak_time,
            "cow_lead_hours": (cow_peak_time - reg_peak_time).total_seconds() / 3600.0,
            "cow_at_reg_crest_cfs": cow_at_reg,
            "lag_factor": cow_at_reg / cow_peak if cow_peak > 0 else np.nan,
            "qual_at_cow_peak": code_at(qual, cow_peak_time),
            "qual_at_reg_crest": code_at(qual, reg_peak_time),
            "censored_in_window": censored,
        })
    table = pd.DataFrame(rows)
    table["bin"] = pd.cut(table["cas_unreg_peak_cfs"], BINS, labels=LABELS)
    return table


def stats_block(table):
    ok = table.dropna(subset=["lag_factor"])
    print("\n" + "=" * 78)
    print("LAG FACTOR BY MAGNITUDE BIN")
    print("=" * 78)
    print("%-8s %4s %8s %8s %8s %8s %8s %8s"
          % ("bin", "n", "median", "mean", "p25", "p75", "min", "max"))
    for b in LABELS:
        s = ok[ok["bin"] == b]["lag_factor"]
        if len(s):
            print("%-8s %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
                  % (b, len(s), s.median(), s.mean(), s.quantile(.25),
                     s.quantile(.75), s.min(), s.max()))
    s = ok["lag_factor"]
    print("%-8s %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
          % ("ALL", len(s), s.median(), s.mean(), s.quantile(.25),
             s.quantile(.75), s.min(), s.max()))

    print("\n" + "=" * 78)
    print("WHAT ACTUALLY DRIVES IT")
    print("=" * 78)
    for label, x, y in [
            ("lag_factor vs event size ", ok["cas_unreg_peak_cfs"], ok["lag_factor"]),
            ("lag_factor vs |lead|     ", ok["cow_lead_hours"].abs(), ok["lag_factor"]),
            ("lead      vs event size  ", ok["cas_unreg_peak_cfs"], ok["cow_lead_hours"].abs())]:
        rho, p = stats.spearmanr(x, y)
        flag = "significant" if p < 0.05 else "NOT significant"
        print("   %s rho %+.3f  p=%.4f   %s" % (label, rho, p, flag))
    print("\n   Lead time explains the lag factor better than size does, and\n"
          "   lead is not itself tied to size. Read the tail value as a\n"
          "   property of these seven storms, not of largeness.")

    print("\n" + "=" * 78)
    print("EVENTS WHOSE COWEEMAN PEAK IS RATING-LIMITED (peak pinned near 3,400)")
    print("=" * 78)
    lost = table[table["lag_factor"].isna() | table["censored_in_window"]]
    if len(lost):
        for _, r in lost.sort_values("cas_unreg_peak_cfs", ascending=False).iterrows():
            print("   %s  CAS unreg %9s  cow peak %7s  code %-3d %s"
                  % (r["event_time"].strftime("%d %b %Y"),
                     format(int(r["cas_unreg_peak_cfs"]), ","),
                     format(int(r["cow_peak_cfs"]), ","),
                     r["qual_at_reg_crest"],
                     CODE_MEANING.get(r["qual_at_reg_crest"], "?")
                     + ("" if np.isfinite(r["lag_factor"])
                        else "  [no lag factor]")))
    print("   -> these are events where the tributary was LARGEST, so every\n"
          "      statistic above is biased LOW by their absence.")


def candidates(table):
    ok = table.dropna(subset=["lag_factor"])
    tail = ok[ok["bin"] == ">60k"]["lag_factor"]
    mid = ok[ok["bin"] == "40-60k"]["lag_factor"]
    opts = [
        ("tail median", tail.median(), "n=%d, the >60k bin" % len(tail)),
        ("tail mean", tail.mean(), "n=%d, pulled up by one 0.81" % len(tail)),
        ("tail p25", tail.quantile(.25), "conservative within the tail"),
        ("tail p75", tail.quantile(.75), "generous within the tail"),
        ("40-60k median", mid.median(), "n=%d, better sampled, less extrapolated" % len(mid)),
        ("all-events median", ok["lag_factor"].median(), "n=%d, ignores the decline" % len(ok)),
    ]
    print("\n" + "=" * 78)
    print("CANDIDATE VALUES AND WHAT EACH COSTS AT AEP=%.4f" % TARGET_AEP)
    print("=" * 78)
    print("%-20s %7s %13s %11s   %s" % ("basis", "value", "local cfs", "combined", "note"))
    base = None
    for name, v, note in opts:
        local = CAS_UNREG_AT_TARGET * LOCAL_RATIO * v
        comb = CAS_REG_AT_TARGET + local
        base = comb if base is None else base
        print("%-20s %7.3f %13s %11s   %s"
              % (name, v, format(int(local), ","), format(int(comb), ","), note))
    lo = CAS_REG_AT_TARGET + CAS_UNREG_AT_TARGET * LOCAL_RATIO * tail.min()
    hi = CAS_REG_AT_TARGET + CAS_UNREG_AT_TARGET * LOCAL_RATIO * tail.max()
    print("\n   full observed tail range %.3f-%.3f spans %s to %s cfs (%.1f%%)"
          % (tail.min(), tail.max(), format(int(lo), ","), format(int(hi), ","),
             100 * (hi - lo) / lo))


def plot(table):
    ok = table.dropna(subset=["lag_factor"])
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax1, ax2, ax3, ax4 = axes.ravel()

    for b in LABELS:
        s = ok[ok["bin"] == b]
        ax1.scatter(s["cas_unreg_peak_cfs"] / 1000.0, s["lag_factor"], s=42,
                    color=C_BIN[b], edgecolor="0.25", lw=0.5, label=b, zorder=3)
    for b in LABELS:
        s = ok[ok["bin"] == b]["lag_factor"]
        lo, hi = {"20-40k": (20, 40), "40-60k": (40, 60),
                  ">60k": (60, ok["cas_unreg_peak_cfs"].max() / 1000.0)}[b]
        ax1.hlines(s.median(), lo, hi, color=C_BIN[b], lw=2.5, ls="--", zorder=4)
    ax1.set_xscale("log")
    ax1.set_xlabel("Castle Rock unregulated peak (1000 cfs)")
    ax1.set_ylabel("lag factor")
    ax1.set_title("Lag factor vs event size\ndashed = bin median")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.scatter(ok["cow_lead_hours"], ok["lag_factor"],
                c=[C_BIN[b] for b in ok["bin"]], s=42, edgecolor="0.25", lw=0.5)
    rho, p = stats.spearmanr(ok["cow_lead_hours"].abs(), ok["lag_factor"])
    ax2.set_xlabel("Coweeman peak lead, hours (negative = leads the crest)")
    ax2.set_ylabel("lag factor")
    ax2.set_title("Lag factor vs lead time\nrho %+.3f, p=%.4f -- the real driver" % (rho, p))
    ax2.grid(alpha=0.3)

    data = [ok[ok["bin"] == b]["lag_factor"].values for b in LABELS]
    try:                      # matplotlib >= 3.9 renamed labels -> tick_labels
        ax3.boxplot(data, tick_labels=LABELS, showmeans=True, widths=0.55)
    except TypeError:
        ax3.boxplot(data, labels=LABELS, showmeans=True, widths=0.55)
    for i, d in enumerate(data, start=1):
        ax3.scatter(np.random.normal(i, 0.05, len(d)), d, s=26, alpha=0.75,
                    color=C_BIN[LABELS[i - 1]], edgecolor="0.25", lw=0.4, zorder=3)
    ax3.set_ylabel("lag factor")
    ax3.set_title("Distribution by bin\nbox = quartiles, triangle = mean, dots = events")
    ax3.grid(alpha=0.3, axis="y")

    tail = np.sort(ok[ok["bin"] == ">60k"]["lag_factor"].values)
    ax4.step(tail, np.arange(1, len(tail) + 1) / len(tail), where="post",
             color="#1a4f8a", lw=2)
    ax4.scatter(tail, np.arange(1, len(tail) + 1) / len(tail), color="#1a4f8a", s=34)
    for v, lab, c in [(np.median(tail), "median %.2f" % np.median(tail), "#b7410e"),
                      (tail.mean(), "mean %.2f" % tail.mean(), "#4c8c4a")]:
        ax4.axvline(v, color=c, ls="--", lw=1.6, label=lab)
    ax4.set_xlim(0, 1)
    ax4.set_xlabel("lag factor")
    ax4.set_ylabel("cumulative fraction")
    ax4.set_title(">60k bin only, all %d events\nthe whole basis for the tail value"
                  % len(tail))
    ax4.grid(alpha=0.3)
    ax4.legend(fontsize=9)

    fig.suptitle("Lag factor evidence -- Coweeman flow at the regulated Castle Rock "
                 "crest, over its own event peak", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


def main():
    print("Reading...")
    cow, qual = read_ecology_cache(CACHE_DIR, label="Coweeman")
    unreg = read_dss_series(DSS_PATH, UNREG_PATHNAME, "Castle Rock unreg :")
    reg = read_dss_series(DSS_PATH, REG_PATHNAME, "Castle Rock reg   :")
    lo = max(cow.index.min(), unreg.index.min())
    hi = min(cow.index.max(), unreg.index.max())
    cow, qual = cow.loc[lo:hi], qual.loc[lo:hi]
    unreg, reg = unreg.loc[lo:hi], reg.loc[lo:hi]

    table = build(cow, qual, unreg, reg)
    table.to_csv(OUT_CSV, index=False)
    print("\nWrote %s (%d events)" % (OUT_CSV, len(table)))
    stats_block(table)
    candidates(table)
    plot(table)


if __name__ == "__main__":
    main()
