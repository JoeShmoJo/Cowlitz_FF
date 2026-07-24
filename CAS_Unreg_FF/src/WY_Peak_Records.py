"""
WY_Peak_Records.py

Extract water-year peak records at Castle Rock from the hourly series:

  UNREG:  //CASTLE ROCK/FLOW-UNREG//1HOUR/CAS+ROUTED-DIFF/  (built by
          Build_Hourly_Holdout_Unreg.py, in MOS_Cleaned.dss)
  REG:    Castle Rock hourly USGS flow (obsData.dss)

For each water year, computes:
  * peak 1-hour value and its timestamp, for UNREG and REG
  * peak 1-day value (max of a 24-hr trailing mean) and the timestamp of
    the window END, for UNREG and REG
  * REG-minus-UNREG peak differences (1-hr and 1-day)
  * peak timing offset (REG peak time minus UNREG peak time, hours) so
    we can track how far apart the two peaks land -- large offsets mean
    the pair is probably not the same storm and the WY should be viewed
    with suspicion before using it in the regression.
  * Oct-Mar valid-hour coverage fractions for both series, as a data
    quality screen.

The REG peak here is deliberately the 1-hr max of the USGS hourly
record, NOT the USGS instantaneous peak-flow record, so that REG and
UNREG peaks are computed identically for the regression in
PeakDiff_Storage_Regression.py.

Output: ../output/wy_peak_records.csv. One row per WY where the UNREG record has any
valid data; REG-only WYs (candidates for the regression correction) are
also included, with the UNREG columns blank.
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
unreg_dir = os.path.join(PROJECT_DIR, "output")
output_dir = os.path.join(PROJECT_DIR, "output")

sys.path.insert(0, REPO_ROOT)
UTILS_DIR = os.path.join(PROJECT_DIR, "src", "Cowlitz_Unreg", "Cowlitz")
sys.path.insert(0, UTILS_DIR)
from utilsDSS import HecDss  # noqa: E402

DSS_OBS = os.path.join(root_dir, "obsData.dss")        # regulated hourly
DSS_UNREG = os.path.join(unreg_dir, "MOS_Cleaned.dss")  # unreg hourly

PATH_CAS_REG = "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1Hour/USGS/"
PATH_CAS_UNREG = "//CASTLE ROCK/FLOW-UNREG//1HOUR/CAS+ROUTED-DIFF/"

OUT_CSV = os.path.join(output_dir, "wy_peak_records.csv")
diag_dir = os.path.join(PROJECT_DIR, "diagnostics")
OUT_GAPS_CSV = os.path.join(diag_dir, "wy_missing_windows.csv")

SEASON_MONTHS = [10, 11, 12, 1, 2, 3]  # Oct-Mar, matches the holdout season

# --- missing-window reporting -------------------------------------------
# Every contiguous Oct-Mar gap >= MIN_GAP_HRS in either series is listed
# in OUT_GAPS_CSV with its distance to that WY's computed peak, so the
# significance of each window can be judged. NOTHING is omitted
# automatically -- review the gap table, then populate the two lists
# below to omit manually and re-run.
MIN_GAP_HRS = 3

# Date ranges to MASK (treated as missing) before peaks are computed.
# Use when a suspect window (bad or gap-adjacent data) should not be
# allowed to set a peak, without dropping the whole WY. Each entry:
# ("YYYY-MM-DD HH:MM", "YYYY-MM-DD HH:MM", series) with series one of
# "reg", "unreg", "both". Endpoints inclusive.
EXCLUDE_RANGES = [
    # ("1996-02-05 00:00", "1996-02-07 12:00", "unreg"),
]

# Water years to drop entirely from the output table (peaks judged
# unusable after reviewing the gap report). These WYs then fall through
# to the dS_2day regression fill in Write_SSP_Record.py if a regulated
# peak exists elsewhere.
OMIT_WYS = [
    # 1997,
]
ONE_DAY_HOURS = 24                     # trailing-mean window for 1-day peak
MIN_ONE_DAY_VALID = 20                 # need >= this many valid hrs in the 24-hr
                                       # window for the mean to count

###############################################################################
# FUNCTION DEFINITIONS


def read_series(dss_file, pathname):
    """Read one hourly series from DSS with sentinel handling; returns a
    pandas Series indexed by datetime (may be empty)."""
    dss = HecDss.open(dss_file)
    try:
        df = dss.readDF(pathname)
    finally:
        dss.close()
    if df.empty:
        return pd.Series(dtype=float)
    s = df["value"]
    s.index = pd.to_datetime(s.index)  # utilsDSS may return a plain Index
    s = s.mask((s <= -900.0))
    s = s[~s.index.duplicated(keep="last")].sort_index()
    return s


def water_year(index):
    """Water year (Oct 1 - Sep 30) for a DatetimeIndex."""
    return index.year + (index.month >= 10).astype(int)


def season_only(s):
    """Restrict a series to the Oct-Mar flood season."""
    return s[s.index.month.isin(SEASON_MONTHS)]


def wy_hourly_peak(s):
    """Per-WY 1-hr max and its timestamp. Returns DataFrame indexed by WY
    with columns [peak, peak_time]."""
    s = season_only(s.dropna())
    if s.empty:
        return pd.DataFrame(columns=["peak", "peak_time"])
    wy = water_year(s.index)
    out = {}
    for y in np.unique(wy):
        sub = s[wy == y]
        t = sub.idxmax()
        out[y] = (sub.loc[t], t)
    df = pd.DataFrame.from_dict(out, orient="index",
                                columns=["peak", "peak_time"])
    df.index.name = "WY"
    return df


def wy_oneday_peak(s):
    """Per-WY max of the 24-hr trailing mean (timestamp = window end).
    Windows with fewer than MIN_ONE_DAY_VALID valid hours are ignored."""
    s = s.copy()
    # regular hourly grid so rolling() window == real hours
    s = s.resample("1h").mean()
    roll = s.rolling(ONE_DAY_HOURS, min_periods=MIN_ONE_DAY_VALID).mean()
    roll = season_only(roll.dropna())
    if roll.empty:
        return pd.DataFrame(columns=["peak", "peak_time"])
    wy = water_year(roll.index)
    out = {}
    for y in np.unique(wy):
        sub = roll[wy == y]
        t = sub.idxmax()
        out[y] = (sub.loc[t], t)
    df = pd.DataFrame.from_dict(out, orient="index",
                                columns=["peak", "peak_time"])
    df.index.name = "WY"
    return df


def find_gaps(s, label):
    """List contiguous Oct-Mar missing windows >= MIN_GAP_HRS.
    Returns a DataFrame [WY, series, gap_start, gap_end, gap_hrs]."""
    if s.dropna().empty:
        return pd.DataFrame(
            columns=["WY", "series", "gap_start", "gap_end", "gap_hrs"])
    idx = pd.date_range(s.index.min(), s.index.max(), freq="1h")
    full = season_only(s.reindex(idx))
    miss = full.isna()
    # group consecutive hourly timestamps (season gaps between Mar->Oct
    # are excluded because season_only removed those rows entirely)
    grp = (
        (~miss) | (miss.index.to_series().diff() > pd.Timedelta(hours=1))
    ).cumsum()
    rows = []
    for _, block in full[miss].groupby(grp[miss]):
        t0, t1 = block.index[0], block.index[-1]
        hrs = (t1 - t0).total_seconds() / 3600.0 + 1
        if hrs >= MIN_GAP_HRS:
            rows.append({"WY": int(water_year(block.index)[0]),
                         "series": label, "gap_start": t0, "gap_end": t1,
                         "gap_hrs": hrs})
    return pd.DataFrame(rows)


def apply_exclusions(reg, unreg):
    """Mask EXCLUDE_RANGES in the requested series (logged)."""
    for start, end, which in EXCLUDE_RANGES:
        t0, t1 = pd.Timestamp(start), pd.Timestamp(end)
        for label, s in (("reg", reg), ("unreg", unreg)):
            if which in (label, "both"):
                n = int(s.loc[t0:t1].notna().sum())
                s.loc[t0:t1] = np.nan
                print(f"  EXCLUDED [{label}] {start} -> {end}: "
                      f"{n} values masked")
    return reg, unreg


def wy_coverage(s):
    """Fraction of Oct-Mar hours with valid data, per WY."""
    idx = pd.date_range(s.index.min(), s.index.max(), freq="1h")
    full = s.reindex(idx)
    full = season_only(full)
    wy = water_year(full.index)
    grp = full.groupby(wy)
    return (grp.count() / grp.size()).rename("coverage")


def build_table(reg, unreg):
    """Assemble the per-WY peak table."""
    reg_1h = wy_hourly_peak(reg)
    unreg_1h = wy_hourly_peak(unreg)
    reg_1d = wy_oneday_peak(reg)
    unreg_1d = wy_oneday_peak(unreg)
    cov_reg = wy_coverage(reg)
    cov_unreg = wy_coverage(unreg) if not unreg.dropna().empty else pd.Series(dtype=float)

    years = sorted(set(reg_1h.index) | set(unreg_1h.index))
    rows = []
    for y in years:
        r = {}
        r["WY"] = y
        r["reg_peak_1hr"] = reg_1h["peak"].get(y, np.nan)
        r["reg_peak_1hr_time"] = reg_1h["peak_time"].get(y, pd.NaT)
        r["unreg_peak_1hr"] = unreg_1h["peak"].get(y, np.nan)
        r["unreg_peak_1hr_time"] = unreg_1h["peak_time"].get(y, pd.NaT)
        r["reg_peak_1day"] = reg_1d["peak"].get(y, np.nan)
        r["reg_peak_1day_time"] = reg_1d["peak_time"].get(y, pd.NaT)
        r["unreg_peak_1day"] = unreg_1d["peak"].get(y, np.nan)
        r["unreg_peak_1day_time"] = unreg_1d["peak_time"].get(y, pd.NaT)

        r["reg_minus_unreg_1hr"] = r["reg_peak_1hr"] - r["unreg_peak_1hr"]
        r["reg_minus_unreg_1day"] = r["reg_peak_1day"] - r["unreg_peak_1day"]

        if pd.notna(r["reg_peak_1hr_time"]) and pd.notna(r["unreg_peak_1hr_time"]):
            dt = r["reg_peak_1hr_time"] - r["unreg_peak_1hr_time"]
            r["peak_offset_hrs"] = dt.total_seconds() / 3600.0
        else:
            r["peak_offset_hrs"] = np.nan

        r["reg_coverage_octmar"] = round(float(cov_reg.get(y, np.nan)), 3)
        r["unreg_coverage_octmar"] = round(float(cov_unreg.get(y, np.nan)), 3) \
            if len(cov_unreg) else np.nan
        rows.append(r)
    return pd.DataFrame(rows).set_index("WY")


###############################################################################
# MAIN


def main():
    print(f"Reading regulated hourly from {DSS_OBS}")
    reg = read_series(DSS_OBS, PATH_CAS_REG)
    print(f"  {len(reg)} values")
    print(f"Reading unreg hourly from {DSS_UNREG}")
    unreg = read_series(DSS_UNREG, PATH_CAS_UNREG)
    print(f"  {len(unreg)} values")

    if EXCLUDE_RANGES:
        print("Applying manual exclusion ranges...")
        reg, unreg = apply_exclusions(reg, unreg)

    # --- missing-window report (identification only; nothing omitted) ---
    gaps = pd.concat([find_gaps(reg, "reg"), find_gaps(unreg, "unreg")],
                     ignore_index=True)

    table = build_table(reg, unreg)

    if len(gaps):
        # distance from each gap to that WY's computed 1-hr peak, per series
        def _dist(row):
            col = f"{row['series']}_peak_1hr_time"
            t = table[col].get(row["WY"], pd.NaT)
            if pd.isna(t):
                return np.nan
            return round(min(abs((t - row["gap_start"]).total_seconds()),
                             abs((t - row["gap_end"]).total_seconds()))
                         / 3600.0, 1)
        gaps["hrs_gap_to_peak"] = gaps.apply(_dist, axis=1)
        gaps = gaps.sort_values(["WY", "series", "gap_start"])
    gaps.to_csv(OUT_GAPS_CSV, index=False)
    print(f"Missing-window report: {len(gaps)} gaps >= {MIN_GAP_HRS} hrs "
          f"-> {OUT_GAPS_CSV}")

    # per-WY gap summary columns for at-a-glance significance checks
    for label in ("reg", "unreg"):
        sub = gaps[gaps["series"] == label] if len(gaps) else gaps
        if len(sub):
            g = sub.groupby("WY")["gap_hrs"]
            table[f"{label}_n_gaps"] = g.count().reindex(table.index)
            table[f"{label}_max_gap_hrs"] = g.max().reindex(table.index)
            near = sub.groupby("WY")["hrs_gap_to_peak"].min()
            table[f"{label}_nearest_gap_to_peak_hrs"] = \
                near.reindex(table.index)
        else:
            table[f"{label}_n_gaps"] = 0
            table[f"{label}_max_gap_hrs"] = np.nan
            table[f"{label}_nearest_gap_to_peak_hrs"] = np.nan
    table[["reg_n_gaps", "unreg_n_gaps"]] = \
        table[["reg_n_gaps", "unreg_n_gaps"]].fillna(0).astype(int)

    if OMIT_WYS:
        dropped = [y for y in OMIT_WYS if y in table.index]
        table = table.drop(index=dropped)
        print(f"OMITTED WYs (manual): {dropped}")

    table.to_csv(OUT_CSV)
    print(f"Wrote {len(table)} WY rows -> {OUT_CSV}")

    both = table.dropna(subset=["reg_minus_unreg_1hr"])
    print(f"\nWYs with both reg and unreg 1-hr peaks: {len(both)}")
    if len(both):
        off = both["peak_offset_hrs"].abs()
        print(f"  peak timing offset |hrs|: median {off.median():.0f}, "
              f"max {off.max():.0f}")
        flagged = both[off > 72]
        if len(flagged):
            print(f"  WARNING: {len(flagged)} WY(s) with peaks > 72 hrs apart "
                  f"(likely different storms): {list(flagged.index)}")


if __name__ == "__main__":
    main()
