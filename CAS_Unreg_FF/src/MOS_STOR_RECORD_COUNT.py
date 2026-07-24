"""
MOS_Daily_Record_Count.py

QC helper for the Castle Rock unreg workflow: counts how many valid
MOS/STOR values are present in each calendar day and writes that count
as a new daily DSS record. Useful for spotting days with partial or
missing telemetry before the despike/smooth step in
Clean_MOS_Holdout.py runs.

What this script does:

1. Reads MOS STOR from the input DSS file (same sentinel handling as
   Clean_MOS_Holdout.py: values <= -900 are treated as missing).
2. Groups the hourly record by calendar day and counts how many
   non-missing values fall in each day (0-24 for an hourly record; a
   full day should be 24, anything less flags a gap).
3. Reindexes to a continuous daily grid over the record's date range so
   days with NO values at all in the DSS record (not just missing/
   sentinel values, but no rows whatsoever) are written as 0 rather
   than skipped.
4. Writes the result back to obsData.dss as a new daily pathname, units
   "Count".
"""

import os
import sys

import pandas as pd

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path
USE_REFERENCE_DATA = True                            # True = ref_data sample run

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")

if USE_REFERENCE_DATA:
    root_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_in")
    output_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_out")
else:
    root_dir = os.path.join(PROJECT_DIR, "data")
    output_dir = os.path.join(PROJECT_DIR, "output")

sys.path.insert(0, REPO_ROOT)
UTILS_DIR = os.path.join(PROJECT_DIR, "src", "Cowlitz_Unreg", "Cowlitz")
sys.path.insert(0, UTILS_DIR)
from utilsDSS import HecDss  # noqa: E402  (project DSS wrapper, handles gaps)

DSS_IN = os.path.join(root_dir, "obsData.dss")
# Convention: obsData.dss holds source records only (observed + hand-cleaned);
# everything a script writes -- including this QC count -- goes to the output
# DSS. Older runs wrote the count back into obsData; that copy can be deleted.
DSS_OUT = os.path.join(output_dir, "MOS_Cleaned.dss")

PATH_MOS_STOR = "//MOS/STOR//1HOUR/CWMS/"
PATH_OUT_COUNT = "//MOS/STOR-COUNT//1DAY/CWMS/"

SENTINEL_THRESH = -900.0  # values <= this are treated as missing, matches
                          # the sentinel handling in Clean_MOS_Holdout.py

###############################################################################
# FUNCTION DEFINITIONS


def load_stor():
    """Read MOS STOR from DSS_IN, mask sentinel/missing values, dedup+sort."""
    dss = HecDss.open(DSS_IN)
    try:
        df = dss.readDF(PATH_MOS_STOR)
        if df.empty:
            raise RuntimeError(f"No data for {PATH_MOS_STOR} in {DSS_IN}")
        s = df["value"]
        s = s.mask(s <= SENTINEL_THRESH)
        s = s[~s.index.duplicated(keep="last")].sort_index()
        print(f"  read STOR: {len(s)} values, {s.index[0]} -> {s.index[-1]}, "
              f"{int(s.isna().sum())} missing")
    finally:
        dss.close()
    return s


def daily_record_count(series):
    """
    Count non-missing values per calendar day, reindexed to a continuous
    daily grid so days with no rows at all in the DSS record show 0
    rather than being dropped.
    """
    counts = series.notna().groupby(series.index.normalize()).sum()
    full_index = pd.date_range(counts.index.min(), counts.index.max(), freq="D")
    counts = counts.reindex(full_index, fill_value=0).astype(float)
    print(f"  {len(counts)} days, {int((counts < 24).sum())} days with "
          f"fewer than 24 values, {int((counts == 0).sum())} days with none")
    return counts


def write_count(counts):
    """
    Write the daily count series to DSS_OUT as PATH_OUT_COUNT, units Count.
    Uses type "PER-AVER" -- utilsDSS's VALID_TYPE_STR doesn't include
    "PER-CUM", and PER-AVER is the type Clean_MOS_Holdout.py already uses
    for its other 1DAY writes (quick_unreg_daily, local_daily).
    """
    dss = HecDss.open(DSS_OUT)
    try:
        ok = dss.writeSeries(counts, PATH_OUT_COUNT, "Count", "PER-AVER")
        print(f"    wrote {PATH_OUT_COUNT}  ({len(counts)} values)"
              if ok else f"    FAILED {PATH_OUT_COUNT}")
    finally:
        dss.close()


def main():
    print("Loading MOS STOR...")
    stor = load_stor()

    print("Counting records per day...")
    counts = daily_record_count(stor)

    print("Writing daily count record to DSS...")
    write_count(counts)

    print("Done.")


###############################################################################
# MAIN

if __name__ == "__main__":
    main()