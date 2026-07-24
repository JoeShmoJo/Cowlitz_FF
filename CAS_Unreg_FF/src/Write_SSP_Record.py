"""
Write_SSP_Record.py

Assemble the final unregulated annual-maximum record at Castle Rock and
write it to a DSS file for import into HEC-SSP. Run this LAST, after:
    1. Build_Hourly_Holdout_Unreg.py   (hourly unreg -> MOS_Cleaned.dss)
    2. WY_Peak_Records.py              (wy_peak_records.csv)
    3. PeakDiff_Storage_Regression.py  (unreg_peak_estimates.csv)
    4. Unreg_Durations_MassBalance.py  (unreg_durations_massbalance.csv)

Sources per duration:
    Peak      unreg_peak_1hr from the hourly record (holdout WYs);
              gap WYs filled from the adopted dS_2day regression
              estimates (unreg_peak_estimates.csv). Source column tags
              which is which.
    One_day   unreg_peak_1day from the hourly record (holdout WYs only;
              no regression fill is defined for 1-day).
    Three_Day daily mass balance (CAS daily + daily MOS holdout from
    Five_Day  the CWMS-CLEAN daily means), restricted to WYs whose
              flood season is complete enough (MAX_SEASON_MISSING_DAYS).

Outputs:
    ../output/wy_record_ssp.csv    audit table (one row per WY, all
                                   durations + source tags)
    ../output/CAS_Unreg_SSP.dss    IR-YEAR records for SSP, one per
                                   duration, F part <DUR>-MOS-HOLDOUT

DSS conventions follow the proven Write_SSP_Inputs.py pattern: 6-part
pathname with the duration folded into the F part (a 7-part template
makes pydsstools silently drop the tag), interval -1, values stamped
Sep 30 of the water year, INST-VAL, CFS.
"""

import os

import numpy as np
import pandas as pd
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
output_dir = os.path.join(PROJECT_DIR, "output")

PEAKS_CSV = os.path.join(output_dir, "wy_peak_records.csv")
ESTIMATES_CSV = os.path.join(output_dir, "unreg_peak_estimates.csv")
MASSBAL_CSV = os.path.join(output_dir, "unreg_durations_massbalance.csv")

OUT_CSV = os.path.join(output_dir, "wy_record_ssp.csv")
OUT_DSS = os.path.join(output_dir, "CAS_Unreg_SSP.dss")

# WYs from the mass-balance table qualify for 3/5-day only if their
# flood season has at most this many missing days.
MAX_SEASON_MISSING_DAYS = 0

# Screen for the hourly-derived peaks: require at least this Oct-Mar
# coverage fraction in the unreg record (column from WY_Peak_Records).
MIN_UNREG_COVERAGE = 0.9

# F part tag identifying the methodology; duration is prefixed per the
# 6-part-pathname rule.
DUR_F = {"Peak": "PEAK", "One_day": "1-DAY",
         "Three_Day": "3-DAY", "Five_Day": "5-DAY"}
PATH_TMPL = "//CASTLEROCK/FLOW-UNREG//IR-YEAR/{durf}-MOS-HOLDOUT/"

###############################################################################
# FUNCTION DEFINITIONS


def assemble_record(peaks, estimates, massbal):
    """Merge the three sources into one WY table with source tags."""
    years = sorted(set(peaks.index) | set(estimates.index) | set(massbal.index))
    rows = []
    for wy in years:
        r = {"WY": wy, "Peak": np.nan, "Peak_Source": "",
             "One_day": np.nan, "One_day_Source": "",
             "Three_Day": np.nan, "Five_Day": np.nan,
             "Durations_Source": ""}

        # --- Peak and One_day from the hourly record ---
        if wy in peaks.index:
            p = peaks.loc[wy]
            cov = p.get("unreg_coverage_octmar", np.nan)
            cov_ok = np.isfinite(cov) and cov >= MIN_UNREG_COVERAGE
            if np.isfinite(p.get("unreg_peak_1hr", np.nan)) and cov_ok:
                r["Peak"] = p["unreg_peak_1hr"]
                r["Peak_Source"] = "hourly_holdout"
            if np.isfinite(p.get("unreg_peak_1day", np.nan)) and cov_ok:
                r["One_day"] = p["unreg_peak_1day"]
                r["One_day_Source"] = "hourly_holdout"

        # --- Peak gap fill from the adopted regression ---
        if not np.isfinite(r["Peak"]) and wy in estimates.index:
            e = estimates.loc[wy]
            if np.isfinite(e.get("unreg_peak_est", np.nan)):
                r["Peak"] = e["unreg_peak_est"]
                r["Peak_Source"] = "dS2day_regression"

        # --- 3/5-day from daily mass balance ---
        if wy in massbal.index:
            m = massbal.loc[wy]
            season_ok = m.get("flood_season_missing", np.inf) \
                <= MAX_SEASON_MISSING_DAYS
            if season_ok:
                for c in ("Three_Day", "Five_Day"):
                    if np.isfinite(m.get(c, np.nan)):
                        r[c] = m[c]
                        r["Durations_Source"] = "daily_massbalance"
        rows.append(r)
    df = pd.DataFrame(rows).set_index("WY")
    return df[df[["Peak", "One_day", "Three_Day", "Five_Day"]]
              .notna().any(axis=1)]


def write_dss(table):
    """Write one IR-YEAR record per duration, Sep 30 stamped."""
    with HecDss.Open(OUT_DSS, version=6) as dss:
        for col, durf in DUR_F.items():
            s = table[col].dropna()
            s = s[np.isfinite(s)]
            if s.empty:
                print(f"  {col}: no valid values; skipping")
                continue
            path = PATH_TMPL.format(durf=durf)
            assert path.count("/") == 7, (
                f"Malformed DSS pathname (needs exactly 6 parts): {path}")
            tsc = TimeSeriesContainer()
            tsc.pathname = path
            tsc.interval = -1
            tsc.times = [pd.Timestamp(f"{int(wy)}-09-30").to_pydatetime()
                         for wy in s.index]
            tsc.values = [float(v) for v in s.values]
            tsc.numberValues = len(s)
            tsc.units = "CFS"
            tsc.type = "INST-VAL"
            try:
                dss.deletePathname(path)
            except Exception:
                pass
            dss.put_ts(tsc)
            print(f"  {col}: {len(s)} WYs ({int(s.index.min())}-"
                  f"{int(s.index.max())}) -> {path}")


###############################################################################
# MAIN


def main():
    print(f"Reading {PEAKS_CSV}")
    peaks = pd.read_csv(PEAKS_CSV, index_col="WY")
    print(f"Reading {ESTIMATES_CSV}")
    estimates = pd.read_csv(ESTIMATES_CSV, index_col="WY") \
        if os.path.exists(ESTIMATES_CSV) else pd.DataFrame()
    print(f"Reading {MASSBAL_CSV}")
    massbal = pd.read_csv(MASSBAL_CSV, index_col="WY")

    table = assemble_record(peaks, estimates, massbal)
    table.to_csv(OUT_CSV)
    print(f"\nAssembled {len(table)} WYs -> {OUT_CSV}")
    for col in ("Peak", "One_day", "Three_Day", "Five_Day"):
        n = int(table[col].notna().sum())
        print(f"  {col}: {n} WYs")
    n_reg = int((table["Peak_Source"] == "dS2day_regression").sum())
    print(f"  ({n_reg} peaks from the dS_2day regression fill)")

    print(f"\nWriting {OUT_DSS}")
    write_dss(table)
    print("\nDone. In HEC-SSP, import the *-MOS-HOLDOUT pathnames.")


if __name__ == "__main__":
    main()
