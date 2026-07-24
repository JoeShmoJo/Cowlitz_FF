"""
Build_Hourly_Holdout_Unreg.py (formerly Clean_MOS_Holdout.py)

Castle Rock unreg workflow, simplified version. MOS ELEV is now hand-
cleaned in DSSVue upstream of this script, so the despike/rolling-
median/Savitzky-Golay/event-detection machinery from the earlier
version is gone. What's left:

1. Reads the hand-cleaned MOS ELEV, raw MOS STOR (used only for a data-
   quality count, see step 3), Mayfield outflow, and Castle Rock flow.
2. Converts ELEV to STOR via the official 2014 rating and differences
   it to get the raw hourly holdout (cfs).
3. Processes the raw holdout:
     a. 3-hour centered rolling average.
     b. Trims to Oct-Mar only (rest of the year set to NaN).
     c. Blanks any calendar day where the raw MOS STOR record has fewer
        than MIN_STOR_VALUES_PER_DAY valid (non-sentinel) values --
        i.e. don't trust a holdout computed from a day with sparse
        storage telemetry (including a day with zero STOR rows at all),
        regardless of what ELEV looked like. STOR_COUNT_OVERRIDES can
        force specific date ranges in or out of this check manually.
4. Trims Mayfield and Castle Rock to the processed holdout's valid
   date range -- the holdout is the limiting record.
5. Routes (Mayfield + processed holdout) to Castle Rock, and routes
   Mayfield alone to Castle Rock, but ONLY over stretches where the
   holdout, Mayfield, and Castle Rock are all valid for at least
   MIN_RUN_DAYS continuous days (see build_castle_unreg). Both routing
   inputs are capped so negative values become 100 cfs (the routing
   routine can't take negative flow). The difference between the two
   routed hydrographs is the routed effect of the holdout; that
   difference is added to the observed Castle Rock record to get the
   unregulated estimate.
6. Writes the processed holdout and all routed/derived series to the
   output DSS file.

Routing (Table 3 calibration, lower Cowlitz River model): the holdout
and the Mayfield-alone flow are each routed through three chained
SSARR reaches spanning Mayfield -> Castle Rock along the Cowlitz
mainstem. The Toutle River's own reach ("Tower to Cowlitz+Toutle") is
excluded -- only the Cowlitz mainstem signal is being routed here, not
the Toutle's inflow.
    1. Mayfield_OUT -> Cowlitz R above Toutle R   KTS=5, n=0.1, 5 phases
    2. Cowlitz R above Toutle R -> Cowlitz+Toutle KTS=1, n=0.2, 1 phase
    3. Cowlitz+Toutle -> Castle Rock              KTS=1, n=0.2, 5 phases
"""

import os
import sys

import numpy as np
import pandas as pd

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
root_dir = os.path.join(PROJECT_DIR, "data")
output_dir = os.path.join(PROJECT_DIR, "output")

sys.path.insert(0, REPO_ROOT)
UTILS_DIR = os.path.join(PROJECT_DIR, "src", "Cowlitz_Unreg", "Cowlitz")
sys.path.insert(0, UTILS_DIR)
from utilsDSS import HecDss  # noqa: E402  (project DSS wrapper, handles gaps)
from HydrologicRouting import SsarrReach  # noqa: E402

DSS_IN = os.path.join(root_dir, "obsData.dss")
DSS_OUT = os.path.join(output_dir, "MOS_Cleaned.dss")

# <-- update the F-part if the hand-edited ELEV was saved under a
# different label in DSSVue than the original "CWMS-CLEAN"
PATH_MOS_ELEV_CLEAN = "//MOS/ELEV//1HOUR/CWMS-CLEAN/"
PATH_MOS_STOR = "//MOS/STOR//1HOUR/CWMS/"  # used only for the daily QC count
PATH_MAY_FLOW = "/COWLITZ RIVER BELOW MAYFIELD DAM, WA/14238000/FLOW//1Hour/USGS/"
PATH_CAS_FLOW = "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1Hour/USGS/"

# --- holdout processing ---
CFS_PER_AF_PER_HR = 43560.0 / 3600.0   # 12.1 cfs per (ac-ft/hr)
ROLL_WINDOW_HRS = 3                    # centered rolling average window
SEASON_MONTHS = [10, 11, 12, 1, 2, 3]  # Oct-Mar; rest of year set to NaN
MIN_STOR_VALUES_PER_DAY = 4            # days with fewer valid STOR readings are blanked

# Manual overrides for specific date ranges, applied to the computed
# daily STOR count before the MIN_STOR_VALUES_PER_DAY check -- e.g. a
# day known to be reliable despite a raw telemetry gap (set a high
# count so it passes), or the reverse (set 0 so it's excluded even
# though enough raw values happened to be present). Each entry is
# (start_date, end_date, override_count), both dates inclusive.
STOR_COUNT_OVERRIDES = [
    ("1997-10-01", "1998-10-01", 0), # bad data.
    ("2006-10-01", "2007-10-01", 0), # misses peak, but could be used for reg/unreg relationship
    ("2011-01-10", "2011-01-31", 24),# filled with CDB data
]

# --- routing ---
NEG_FLOW_FLOOR_CFS = 100.0  # routing can't take negative flow; negatives -> this value
ROUTING_TIMESTEP_HRS = 1.0  # hourly data
MIN_RUN_DAYS = 10           # only route stretches with >= this many continuous
                            # valid days (holdout + Mayfield + Castle Rock all
                            # present) -- shorter runs don't give SsarrReach
                            # enough time to settle away from its initial
                            # condition before the segment ends

# Table 3 calibration, lower Cowlitz River model. Only the Cowlitz mainstem
# chain Mayfield -> Castle Rock is used; the Toutle River's own reach
# ("Tower to Cowlitz+Toutle") is excluded since only the holdout signal is
# being routed, not the Toutle's inflow. Each tuple is (kts, n, numSubreaches/phases).
COWLITZ_REACH_PARAMS = [
    (5, 0.1, 5),  # Mayfield_OUT -> Cowlitz R above Toutle R
    (1, 0.2, 1),  # Cowlitz R above Toutle R -> Cowlitz+Toutle
    (1, 0.2, 5),  # Cowlitz+Toutle -> Castle Rock
]

# --- output pathnames (D part left blank; DSS assigns blocks) ---
MISSING_SENTINEL = -902.0  # written in place of gaps so records stay continuous
PATH_OUT = {
    "holdout_raw":             "//MOS/FLOW-HOLDOUT//1HOUR/RAW/",
    "holdout_processed":       "//MOS/FLOW-HOLDOUT//1HOUR/PROCESSED/",
    "routed_may_plus_holdout": "//CASTLE ROCK/FLOW-ROUTED//1HOUR/MAY+HOLDOUT/",
    "routed_may_only":         "//CASTLE ROCK/FLOW-ROUTED//1HOUR/MAY-ONLY/",
    "routed_diff":             "//CASTLE ROCK/FLOW-ROUTED-DIFF//1HOUR/MAYHOLDOUT-MINUS-MAYONLY/",
    "unreg_castle":            "//CASTLE ROCK/FLOW-UNREG//1HOUR/CAS+ROUTED-DIFF/",
}

# --- official MOS elev-storage rating (2014), 1-ft steps 622-780 ft ---
# Implied incremental areas: ~5,750 ac at 630 ft -> ~11,500 ac at 778 ft;
# ~1.30M ac-ft at full pool 778.5. Extrapolated on the edge slope outside
# 622-780 (only matters for surcharge above 780).
OFFICIAL_RATING_ELEV = np.arange(622.0, 781.0, 1.0)

OFFICIAL_RATING_STOR = np.array([
      2579,   7736,  12893,  18050,  23405,  28562,  33719,  38876,  44033,
     49785,  55339,  61091,  66645,  72397,  77950,  83702,  89256,  95008, 100562,
    106512, 112661, 118612, 124760, 130711, 136661, 142810, 148760, 154711, 160859,
    167207, 173554, 179901, 186446, 192793, 199140, 205487, 211835, 218380, 224727,
    231471, 238413, 245157, 251901, 258843, 265587, 272529, 279273, 286016, 292959,
    300297, 307636, 314777, 322116, 329454, 336793, 344132, 351471, 358810, 366149,
    373884, 381620, 389355, 397091, 404826, 412562, 420099, 427834, 435570, 443306,
    451438, 459570, 467702, 475636, 483768, 491901, 500033, 507967, 516099, 524231,
    532760, 541091, 549620, 558148, 566677, 575008, 583537, 592066, 600595, 608925,
    617851, 626777, 635504, 644429, 653157, 662082, 671008, 679735, 688661, 697388,
    706710, 715834, 725157, 734281, 743405, 752727, 761851, 770975, 780297, 789421,
    798942, 808661, 818181, 827702, 837223, 846743, 856264, 865983, 875504, 885024,
    894942, 905057, 914975, 924892, 934809, 944925, 954842, 964760, 974677, 984793,
    995107, 1005619, 1016132, 1026446, 1036958, 1047272, 1057785, 1068297,
  1078611, 1089123, 1100032, 1111140, 1122049, 1133156, 1144065, 1155173,
  1166082, 1177189, 1188099, 1199206, 1210710, 1222413, 1233917, 1245619,
  1257322, 1268826, 1280528, 1292231, 1303735, 1315437,
], dtype=float)

###############################################################################
# FUNCTION DEFINITIONS


def load_inputs():
    """Read ELEV (hand-cleaned), STOR (raw, for QC only), Mayfield, and
    Castle Rock from DSS_IN. Returns dict of Series."""
    data = {}
    dss = HecDss.open(DSS_IN)
    try:
        for key, path in [("elev", PATH_MOS_ELEV_CLEAN), ("stor", PATH_MOS_STOR),
                          ("may", PATH_MAY_FLOW), ("cas", PATH_CAS_FLOW)]:
            df = dss.readDF(path)
            if df.empty:
                raise RuntimeError(f"No data for {path} in {DSS_IN}")
            s = df["value"]
            s = s.mask((s <= -900.0) | (s < -9000.0))  # sentinel handling
            s = s[~s.index.duplicated(keep="last")].sort_index()
            data[key] = s
            print(f"  read {key}: {len(s)} values, {s.index[0]} -> "
                  f"{s.index[-1]}, {int(s.isna().sum())} missing")
    finally:
        dss.close()
    return data


def _interp_extrap(x, xg, yg):
    """np.interp with linear edge-slope extrapolation beyond the grid.

    Extrapolating matters because np.interp's default clamping (or NaN)
    would corrupt holdouts exactly at high pool -- i.e., during floods.
    """
    y = np.interp(x, xg, yg)
    lo_slope = (yg[1] - yg[0]) / (xg[1] - xg[0])
    hi_slope = (yg[-1] - yg[-2]) / (xg[-1] - xg[-2])
    y = np.where(x < xg[0], yg[0] + (x - xg[0]) * lo_slope, y)
    y = np.where(x > xg[-1], yg[-1] + (x - xg[-1]) * hi_slope, y)
    return y


def elev_to_stor(elev_series, elev_grid=OFFICIAL_RATING_ELEV, stor_grid=OFFICIAL_RATING_STOR):
    """Convert an ELEV series to STOR (ac-ft) with the official rating."""
    out = _interp_extrap(elev_series.values, elev_grid, stor_grid)
    return pd.Series(out, index=elev_series.index).where(elev_series.notna())


def hourly_holdout(stor_series):
    """Hourly holdout in cfs from an hourly STOR series (centered diff)."""
    ds = stor_series.diff()  # ac-ft per hour
    return ds * CFS_PER_AF_PER_HR


def daily_stor_counts(stor_series):
    """
    Count of valid (non-missing) STOR values per calendar day, over a
    continuous daily grid spanning the record. Days with zero raw STOR
    rows at all (not just all-sentinel/missing rows) get an explicit 0
    rather than being absent from the result -- absent days would
    otherwise silently pass the MIN_STOR_VALUES_PER_DAY check below.
    """
    counts = stor_series.notna().groupby(stor_series.index.normalize()).sum()
    full_days = pd.date_range(counts.index.min(), counts.index.max(), freq="D")
    return counts.reindex(full_days, fill_value=0)


def apply_stor_count_overrides(stor_counts, overrides=STOR_COUNT_OVERRIDES):
    """
    Apply STOR_COUNT_OVERRIDES to a daily STOR count Series -- for each
    (start_date, end_date, override_count) entry, set every day in that
    inclusive range to override_count. Days outside stor_counts' date
    range are skipped with a warning rather than silently ignored.
    """
    stor_counts = stor_counts.copy()
    for start, end, count in overrides:
        idx = pd.date_range(start, end, freq="D").intersection(stor_counts.index)
        if len(idx) == 0:
            print(f"    WARNING: override {start} -> {end} has no matching "
                  f"days in the STOR record -- skipped")
            continue
        stor_counts.loc[idx] = count
        print(f"    override: {start} -> {end} ({len(idx)} days) count set to {count}")
    return stor_counts


def process_holdout(holdout_raw, stor_counts):
    """
    3-hour centered rolling average -> Oct-Mar only -> blank days with
    too few STOR readings. Each step only removes/smooths; it doesn't
    fill genuinely missing periods beyond the rolling window's reach.
    """
    rolled = holdout_raw.rolling(ROLL_WINDOW_HRS, center=True, min_periods=1).mean()

    in_season = rolled.index.month.isin(SEASON_MONTHS)
    rolled = rolled.where(in_season)

    low_count_days = stor_counts[stor_counts < MIN_STOR_VALUES_PER_DAY].index
    day_of_hour = rolled.index.normalize()
    rolled = rolled.mask(day_of_hour.isin(low_count_days))

    n_total = len(rolled)
    n_kept = int(rolled.notna().sum())
    print(f"    processed holdout: {n_kept}/{n_total} hours kept "
          f"({100.0 * n_kept / max(n_total, 1):.1f}%) after season/QC trim")
    return rolled


def trim_to_holdout_extent(holdout_processed, may, cas):
    """
    Reindex holdout/Mayfield/Castle Rock to the processed holdout's
    valid date range (first_valid_index to last_valid_index) -- the
    holdout is the limiting record for this analysis.
    """
    t0 = holdout_processed.first_valid_index()
    t1 = holdout_processed.last_valid_index()
    if t0 is None:
        raise RuntimeError("Processed holdout has no valid values -- "
                           "check SEASON_MONTHS / MIN_STOR_VALUES_PER_DAY.")
    idx = pd.date_range(t0, t1, freq="h")
    print(f"    trimming to holdout extent: {t0} -> {t1} ({len(idx)} hours)")
    return (holdout_processed.reindex(idx), may.reindex(idx), cas.reindex(idx))


def cap_negative_flows(series, floor=NEG_FLOW_FLOOR_CFS):
    """
    Replace negative flows with `floor` cfs -- the routing routine
    can't take negative input. Uses .mask() (not .where()) so NaN gaps
    are left as NaN: NaN < 0 is False, so .mask() never touches them --
    only true negatives are replaced, and missing data stays missing.
    """
    return series.mask(series < 0, floor)


def find_valid_runs(valid_mask, min_hours):
    """
    Identify maximal contiguous True runs in a boolean, regular-hourly
    Series. Returns a list of (start, end) timestamps for runs spanning
    at least min_hours (inclusive of both endpoints).
    """
    run_id = (valid_mask != valid_mask.shift()).cumsum()
    runs = []
    for _, group in valid_mask.groupby(run_id):
        if not bool(group.iloc[0]):
            continue
        span_hours = (group.index[-1] - group.index[0]) / pd.Timedelta(hours=1) + 1
        if span_hours >= min_hours:
            runs.append((group.index[0], group.index[-1]))
    return runs


def build_cowlitz_reaches():
    """Fresh SsarrReach objects for the three-reach Mayfield -> Castle
    Rock chain (COWLITZ_REACH_PARAMS). Built fresh per call so a route_reach()
    call never reuses another call's internal state (subreachOutflows)."""
    reaches = []
    for kts, n, num_subreaches in COWLITZ_REACH_PARAMS:
        reach = SsarrReach(timestepHrs=ROUTING_TIMESTEP_HRS)
        reach.buildWithKTS(numSubreaches=num_subreaches, n=n, kts=kts)
        reaches.append(reach)
    return reaches


def route_reach(flow_series):
    """
    Route an hourly flow series through the three chained SSARR reaches
    spanning Mayfield -> Castle Rock (Table 3 calibration). flow_series
    should already be non-negative (see cap_negative_flows) -- passed
    with allowNegatives=True since that capping is done explicitly
    upstream rather than left to SsarrReach's own floor-to-1-cfs
    behavior. Returns a pandas Series on the same index as flow_series.
    """
    values = flow_series.values.tolist()
    for reach in build_cowlitz_reaches():
        values = reach.routeHydrograph(values, allowNegatives=True)
    return pd.Series(values, index=flow_series.index)


def build_castle_unreg(holdout_processed, may, cas, min_run_days=MIN_RUN_DAYS):
    """
    Routes (Mayfield + holdout) and Mayfield alone only over stretches
    where the processed holdout, Mayfield, and Castle Rock are all
    valid for at least `min_run_days` continuous days. Shorter runs are
    excluded outright -- SsarrReach initializes each subreach's storage
    to the run's first inflow value, so a short run never gets past
    that artificial start before it ends. Each qualifying run is routed
    with its own fresh reach objects (no state carried between runs or
    across the excluded gaps). Everything outside a qualifying run is
    NaN in the output; no gap-filling or interpolation is used anywhere.
    """
    valid = holdout_processed.notna() & may.notna() & cas.notna()
    runs = find_valid_runs(valid, min_hours=min_run_days * 24)

    idx = holdout_processed.index
    routed_combined = pd.Series(np.nan, index=idx)
    routed_may_only = pd.Series(np.nan, index=idx)

    print(f"    {len(runs)} run(s) of >= {min_run_days} continuous days qualify")
    for t0, t1 in runs:
        holdout_run = holdout_processed.loc[t0:t1]
        may_run = may.loc[t0:t1]

        combined_run = cap_negative_flows(may_run + holdout_run)
        may_capped_run = cap_negative_flows(may_run)

        routed_combined.loc[t0:t1] = route_reach(combined_run).values
        routed_may_only.loc[t0:t1] = route_reach(may_capped_run).values
        n_days = int((t1 - t0) / pd.Timedelta(days=1)) + 1
        print(f"      routed {t0} -> {t1} ({n_days} days)")

    routed_diff = routed_combined - routed_may_only
    unreg_castle = cas + routed_diff

    return routed_combined, routed_may_only, routed_diff, unreg_castle


def write_outputs(holdout_raw, holdout_processed, routed_combined,
                  routed_may_only, routed_diff, unreg_castle):
    """Write the processed holdout and (if available) the routed/derived
    series to the output DSS file via utilsDSS.writeSeries."""
    os.makedirs(output_dir, exist_ok=True)
    jobs = [
        ("holdout_raw", holdout_raw, "cfs", "INST-VAL"),
        ("holdout_processed", holdout_processed, "cfs", "INST-VAL"),
        ("routed_may_plus_holdout", routed_combined, "cfs", "INST-VAL"),
        ("routed_may_only", routed_may_only, "cfs", "INST-VAL"),
        ("routed_diff", routed_diff, "cfs", "INST-VAL"),
        ("unreg_castle", unreg_castle, "cfs", "INST-VAL"),
    ]
    dss = HecDss.open(DSS_OUT)
    try:
        for key, series, units, dtype in jobs:
            if series is None or series.dropna().empty:
                print(f"    skip {key}: not available / empty")
                continue
            s = series.loc[series.first_valid_index():series.last_valid_index()]
            s = s.reindex(pd.date_range(s.index[0], s.index[-1], freq="h"))
            n_missing = int(s.isna().sum())
            s = s.fillna(MISSING_SENTINEL)
            ok = dss.writeSeries(s, PATH_OUT[key], units, dtype)
            print(f"    wrote {PATH_OUT[key]}  ({len(s)} values, "
                  f"{n_missing} as {MISSING_SENTINEL:g})"
                  if ok else f"    FAILED {PATH_OUT[key]}")
    finally:
        dss.close()


def main():
    os.makedirs(output_dir, exist_ok=True)

    print("Loading inputs...")
    data = load_inputs()
    elev, stor, may, cas = data["elev"], data["stor"], data["may"], data["cas"]

    print("Converting hand-cleaned ELEV to STOR via official rating...")
    stor_from_elev = elev_to_stor(elev)

    print("Computing raw hourly holdout...")
    holdout_raw = hourly_holdout(stor_from_elev)

    print("Counting daily STOR values for QC...")
    stor_counts = daily_stor_counts(stor)
    stor_counts = apply_stor_count_overrides(stor_counts)

    print("Processing holdout (3-hr avg, Oct-Mar only, drop low-count days)...")
    holdout_p = process_holdout(holdout_raw, stor_counts)

    print("Trimming Mayfield/Castle Rock to holdout extent...")
    holdout_p, may_t, cas_t = trim_to_holdout_extent(holdout_p, may, cas)

    print("Routing (Mayfield+holdout and Mayfield-only to Castle Rock)...")
    routed_combined, routed_may_only, routed_diff, unreg_castle = \
        build_castle_unreg(holdout_p, may_t, cas_t)

    print("Writing DSS outputs...")
    write_outputs(holdout_raw, holdout_p, routed_combined, routed_may_only,
                  routed_diff, unreg_castle)

    print("Done.")


###############################################################################
# MAIN

if __name__ == "__main__":
    main()