#Coweeman_RegPeak_Timing.py
# -*- coding: utf-8 -*-
"""
How much of the Coweeman's peak is left when the REGULATED Castle Rock
peak arrives -- not the unregulated one.

WHY THIS EXISTS
    coweeman_proportion.csv (#Coweeman_Proportion.py) and everything built
    on it so far measures the Coweeman's recession relative to the
    UNREGULATED Castle Rock peak's timing -- by design, per that script's
    own docstring, since a ratio taken against REGULATED flow carries the
    reservoir's operation in its denominator. That reasoning is correct
    for the MAGNITUDE ratio it was built to measure.

    But the downstream-confluence combination cares about a different
    question: how much of the Coweeman's peak survives until the
    REGULATED Castle Rock peak specifically, since that is the moment the
    combined flow downstream of the confluence is actually being assessed
    at. If regulation delays the Castle Rock peak relative to the
    unregulated peak -- which it can, when the reservoir fills and starts
    passing inflow through after the local peak has already passed -- the
    Coweeman will have receded further by the time the REGULATED peak
    shows up than the unregulated-timing ratio suggests. Asked directly
    whether the unregulated-timing analysis actually answers this; it
    doesn't, on its own.

METHOD
    Same 79-ish independent unregulated events #Coweeman_Proportion.py
    uses (same selection function, same parameters, read fresh here so
    this script has no import dependency on that one). For each event:
      1. Unregulated Castle Rock peak time and magnitude (as before).
      2. REGULATED Castle Rock peak, found by searching a window around
         the unregulated peak time on the REGULATED ResSim series
         (//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/ -- the same DSS
         file's regulated companion to Flow-UNREG).
      3. lag_reg_vs_unreg_hours = regulated peak time minus unregulated
         peak time. Positive means regulation DELAYED the peak.
      4. Coweeman's flow at the REGULATED peak's timestamp (same
         resample-and-lookup approach as #Coweeman_Proportion.py), and
         the ratio of that to the Coweeman's own peak in the window --
         directly comparable to ratio_coincident in coweeman_proportion.csv,
         but anchored to the regulated peak instead of the unregulated one.

    SEARCH WINDOW for the regulated peak: +/-72 hours around the
    unregulated peak. Wide enough to catch a multi-day fill-and-pass-
    through delay without reaching into the next independent event
    (events are >=7 days apart by construction). Any event where the
    regulated peak lands within a few hours of a window edge is flagged
    rather than trusted -- it means 72 hours wasn't enough margin for
    that event specifically.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import re
import glob
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pydsstools.heclib.dss import HecDss

import sys
REPO_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "Modules"))
# The Ecology parser lives in /Modules because four scripts once carried
# copy-pasted copies of it and all four shared the same bug: quality code 254
# ("Rating table exceeded, data will not be reported") was parsed as a
# discharge of 254 cfs. See Modules/ecology_io.py.
from ecology_io import (read_ecology_cache, resample_censor_aware,
                        MISSING_CODES, TRUSTED_CODES, CODE_MEANING)

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
CACHE_DIR = r"../data/coweeman"
DSS_PATH = r"../output/ResSim_WCM_RC.dss"
UNREG_PATHNAME = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"
REG_PATHNAME = "//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "coweeman_regpeak_timing.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coweeman_regpeak_timing.png")

RESAMPLE = "1h"
N_EVENTS = 80
EVENT_MIN_SEPARATION_DAYS = 7
MIN_EVENT_CFS = 20000.0
EVENT_WINDOW_HOURS = 48          # Coweeman-own-peak search window (as before)
REG_SEARCH_WINDOW_HOURS = 72     # regulated-peak search window around unreg peak
MIN_WINDOW_COVERAGE = 0.80

MAGNITUDE_BINS = [20000, 40000, 60000, 200000]
BIN_LABELS = ["20-40k", "40-60k", ">60k"]

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


def read_dss_series(path, pathname):
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
    return pd.Series(values, index=index).dropna().sort_index()


def find_events(series, n, min_separation_days, floor):
    working = series.copy()
    events, span = [], pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span) & (working.index <= when + span)] = np.nan
    return sorted(events)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Reading series...")
    cow, qual = read_ecology_cache(CACHE_DIR, label="Coweeman")
    unreg = read_dss_series(DSS_PATH, UNREG_PATHNAME)
    reg = read_dss_series(DSS_PATH, REG_PATHNAME)
    print("   Coweeman   : %d values, %s to %s" % (len(cow), cow.index.min().date(), cow.index.max().date()))
    print("   Unreg CasRk: %d hours, %s to %s" % (len(unreg), unreg.index.min().date(), unreg.index.max().date()))
    print("   Reg   CasRk: %d hours, %s to %s" % (len(reg), reg.index.min().date(), reg.index.max().date()))

    first = max(cow.index.min(), unreg.index.min(), reg.index.min())
    last = min(cow.index.max(), unreg.index.max(), reg.index.max())
    grid = pd.date_range(first.ceil(RESAMPLE), last.floor(RESAMPLE), freq=RESAMPLE)
    cow_h = resample_censor_aware(cow, qual, RESAMPLE).reindex(grid)
    unreg_h = unreg.resample(RESAMPLE).mean().reindex(grid)
    reg_h = reg.resample(RESAMPLE).mean().reindex(grid)
    print("   overlap    : %s to %s (%d hours)" % (first.date(), last.date(), len(grid)))

    events = find_events(unreg_h, N_EVENTS, EVENT_MIN_SEPARATION_DAYS, MIN_EVENT_CFS)
    print("   %d independent unregulated events above %s cfs" % (len(events), format(int(MIN_EVENT_CFS), ",")))

    half_cow = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
    half_reg = pd.Timedelta(hours=REG_SEARCH_WINDOW_HOURS)
    rows = []
    for when, unreg_peak in events:
        cow_window = slice(when - half_cow, when + half_cow)
        a, b = unreg_h.loc[cow_window], cow_h.loc[cow_window]
        need = int(MIN_WINDOW_COVERAGE * len(a))
        if len(a) == 0 or a.notna().sum() < need or b.notna().sum() < need:
            continue
        at_unreg_peak = b.get(when, np.nan)
        if not np.isfinite(at_unreg_peak):
            continue
        cow_peak_cfs = float(b.max())

        reg_window = slice(when - half_reg, when + half_reg)
        c = reg_h.loc[reg_window]
        if c.notna().sum() < int(MIN_WINDOW_COVERAGE * len(c)):
            continue
        reg_peak_time = c.idxmax()
        reg_peak_cfs = float(c.max())
        at_edge = (reg_peak_time <= reg_window.start + pd.Timedelta(hours=3) or
                   reg_peak_time >= reg_window.stop - pd.Timedelta(hours=3))

        at_reg_peak = cow_h.get(reg_peak_time, np.nan)
        if not np.isfinite(at_reg_peak):
            continue

        rows.append({
            "event_time_unreg": when,
            "cas_unreg_peak_cfs": unreg_peak,
            "reg_peak_time": reg_peak_time,
            "cas_reg_peak_cfs": reg_peak_cfs,
            "lag_reg_vs_unreg_hours": (reg_peak_time - when).total_seconds() / 3600.0,
            "reg_peak_near_window_edge": bool(at_edge),
            "cow_peak_cfs": cow_peak_cfs,
            "cow_at_unreg_peak_cfs": float(at_unreg_peak),
            "cow_at_reg_peak_cfs": float(at_reg_peak),
            "ratio_at_unreg_peak": float(at_unreg_peak) / cow_peak_cfs,
            "ratio_at_reg_peak": float(at_reg_peak) / cow_peak_cfs,
        })

    table = pd.DataFrame(rows)
    table["magnitude_bin"] = pd.cut(table["cas_unreg_peak_cfs"], bins=MAGNITUDE_BINS,
                                     labels=BIN_LABELS, right=False)
    table.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV, "(%d events)" % len(table))

    n_edge = int(table["reg_peak_near_window_edge"].sum())
    if n_edge:
        print("WARNING: %d/%d events had their regulated peak within 3h of the "
              "+/-%dh search window edge -- widen REG_SEARCH_WINDOW_HOURS and "
              "rerun before trusting those rows." % (n_edge, len(table), REG_SEARCH_WINDOW_HOURS))

    print()
    print("Lag (regulated peak time minus unregulated peak time), by magnitude bin:")
    for b, g in table.groupby("magnitude_bin", observed=True):
        lag = g["lag_reg_vs_unreg_hours"]
        print("  %-8s n=%2d  mean=%+.1fh  median=%+.1fh  min=%+.1fh  max=%+.1fh"
              % (b, len(g), lag.mean(), lag.median(), lag.min(), lag.max()))

    print()
    print("Coweeman ratio at unregulated peak vs. at regulated peak, by magnitude bin:")
    for b, g in table.groupby("magnitude_bin", observed=True):
        ru, rr = g["ratio_at_unreg_peak"], g["ratio_at_reg_peak"]
        print("  %-8s n=%2d  unreg-timing mean=%.3f median=%.3f  |  reg-timing mean=%.3f median=%.3f"
              % (b, len(g), ru.mean(), ru.median(), rr.mean(), rr.median()))

    tail = table[table["magnitude_bin"] == ">60k"]
    print()
    print(tail[["event_time_unreg", "cas_unreg_peak_cfs", "lag_reg_vs_unreg_hours",
                "reg_peak_near_window_edge", "ratio_at_unreg_peak", "ratio_at_reg_peak"]]
          .sort_values("cas_unreg_peak_cfs", ascending=False).to_string(index=False))

    # -- plot: ratio at unreg-peak vs reg-peak timing, by magnitude --
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.scatter(table["cas_unreg_peak_cfs"], table["ratio_at_unreg_peak"],
               color="#1a4f8a", s=30, alpha=0.7, label="Ratio at UNREGULATED peak timing")
    ax.scatter(table["cas_unreg_peak_cfs"], table["ratio_at_reg_peak"],
               color="#b7410e", s=30, alpha=0.7, marker="^", label="Ratio at REGULATED peak timing")
    for _, row in table.iterrows():
        ax.plot([row["cas_unreg_peak_cfs"]] * 2,
                [row["ratio_at_unreg_peak"], row["ratio_at_reg_peak"]],
                color="gray", lw=0.6, alpha=0.5, zorder=0)
    ax.axhline(0.80, color="green", lw=1.2, ls="--", label="2009-era borrowed figure (0.80)")
    ax.set_xscale("log")
    ax.set_xlabel("Unregulated Castle Rock peak (cfs)")
    ax.set_ylabel("Coweeman flow / Coweeman's own peak")
    ax.set_title("Coweeman recession: at unregulated vs. regulated Castle Rock peak timing")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


if __name__ == "__main__":
    main()
