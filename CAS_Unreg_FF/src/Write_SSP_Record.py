"""
Write_SSP_Record.py

Assemble the final unregulated annual-maximum record at Castle Rock and
write it to a DSS file for import into HEC-SSP. Run this LAST, after:
    1. Build_Hourly_Holdout_Unreg.py   (hourly unreg -> MOS_Cleaned.dss)
    2. WY_Peak_Records.py              (wy_peak_records.csv)
    3. PeakDiff_Storage_Regression.py  (unreg_peak_estimates.csv)
    4. Unreg_Durations_MassBalance.py  (unreg_durations_massbalance.csv)

Source rules (stated 24 Jul 2026):
    Pre-1968 (<= PRE_REG_LAST_WY, before Mossyrock closure):
        Peak            USGS peak flow record (observed instantaneous
                        annual peaks; unregulated pre-dam)
        One/Three/Five  rolling 1/3/5-day maxima computed directly from
                        the USGS daily record, for WYs whose Oct-Mar
                        season is complete enough
    Post-1968:
        Peak      calculated hourly unreg (1-hr max), else the
                  change-in-daily-storage (dS_2day) regression estimate
                  (unreg_peak_estimates.csv)
        One_day   1-day average of the hourly unreg record; if hourly
                  doesn't exist for the WY, the one-day max of the
                  daily unreg record (mass balance One_day)
        Three_Day the unreg daily averages (mass balance), restricted
        Five_Day  to WYs whose flood season is complete enough
                  (MAX_SEASON_MISSING_DAYS)
    Nothing is read from the Cowlitz_FF_DataPrep archive.

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

# WYs from the mass-balance table (and pre-reg WYs from the daily
# record) qualify for durations only if their Oct-Mar flood season has
# at most this many missing days.
MAX_SEASON_MISSING_DAYS = 0

# --- pre-regulation component (WY <= PRE_REG_LAST_WY) ---
PRE_REG_LAST_WY = 1968  # last pre-regulation WY (Mossyrock closure Dec 1968;
                        # WY1927-1968 treated as unregulated, matching the
                        # 2009 study / archived Build_Simplified convention)
# USGS observed instantaneous annual peaks (active project input;
# unregulated by definition for the pre-dam years)
PRE_REG_PEAKS_CSV = os.path.join(PROJECT_DIR, "data",
                                 "CastleRock_USGS_peaks.csv")
# USGS daily record at Castle Rock for the pre-reg 1/3/5-day durations
DAILY_DSS = os.path.join(PROJECT_DIR, "data", "obsData.dss")
PATH_CAS_DAILY = "/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW//1DAY/USGS/"

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


def read_daily():
    """USGS daily flow at Castle Rock as a Series (sentinels dropped)."""
    ts = None
    with HecDss.Open(DAILY_DSS, version=6) as dss:
        ts = dss.read_ts(PATH_CAS_DAILY)
    s = pd.Series(ts.values, index=pd.to_datetime(ts.pytimes))
    s = pd.to_numeric(s, errors="coerce")
    s = s[s > -900]
    s.index = s.index.normalize()
    s = s[~s.index.duplicated(keep="last")].sort_index()
    return s


def prereg_peaks():
    """USGS instantaneous annual peaks for pre-regulation WYs."""
    if not os.path.exists(PRE_REG_PEAKS_CSV):
        print(f"  WARNING: {PRE_REG_PEAKS_CSV} not found; "
              "no pre-reg peaks will be written")
        return pd.Series(dtype=float)
    df = pd.read_csv(PRE_REG_PEAKS_CSV)
    df = df[df["WY"] <= PRE_REG_LAST_WY]
    return df.set_index("WY")["Peak_cfs"].astype(float)


def prereg_durations(daily):
    """1/3/5-day WY maxima from the USGS daily record for pre-reg WYs
    with a complete-enough Oct-Mar season. Returns DataFrame indexed by
    WY with One_day/Three_Day/Five_Day."""
    wy = daily.index.year + (daily.index.month >= 10).astype(int)
    rows = []
    for y in sorted(set(wy)):
        if y > PRE_REG_LAST_WY:
            continue
        grp = daily[wy == y]
        season = pd.date_range(f"{y-1}-10-01", f"{y}-03-31", freq="1D")
        missing = len(season) - grp.index.isin(season).sum()
        if missing > MAX_SEASON_MISSING_DAYS:
            print(f"    pre-reg WY{y}: {missing} flood-season days "
                  "missing -- durations skipped")
            continue
        r = {"WY": y}
        for label, win in (("One_day", 1), ("Three_Day", 3),
                           ("Five_Day", 5)):
            roll = grp.rolling(win, min_periods=win).mean()
            r[label] = roll.max()
        rows.append(r)
    if not rows:
        return pd.DataFrame(
            columns=["One_day", "Three_Day", "Five_Day"]).rename_axis("WY")
    return pd.DataFrame(rows).set_index("WY")


def assemble_record(peaks, estimates, massbal,
                    prereg_pk=None, prereg_dur=None):
    """Merge all sources into one WY table with source tags."""
    if prereg_pk is None:
        prereg_pk = pd.Series(dtype=float)
    if prereg_dur is None:
        prereg_dur = pd.DataFrame(
            columns=["One_day", "Three_Day", "Five_Day"]).rename_axis("WY")
    years = sorted(set(peaks.index) | set(estimates.index)
                   | set(massbal.index) | set(prereg_pk.index)
                   | set(prereg_dur.index))
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

        # --- pre-regulation WYs: observed USGS peaks + daily durations ---
        if wy <= PRE_REG_LAST_WY:
            if wy in prereg_pk.index and np.isfinite(prereg_pk[wy]):
                r["Peak"] = prereg_pk[wy]
                r["Peak_Source"] = "usgs_peak_prereg"
            if wy in prereg_dur.index:
                d = prereg_dur.loc[wy]
                for c in ("One_day", "Three_Day", "Five_Day"):
                    if np.isfinite(d.get(c, np.nan)):
                        r[c] = d[c]
                r["One_day_Source"] = "usgs_daily_prereg" \
                    if np.isfinite(r["One_day"]) else r["One_day_Source"]
                r["Durations_Source"] = "usgs_daily_prereg"
            rows.append(r)
            continue

        # --- 1/3/5-day from the daily unreg (mass balance) ---
        if wy in massbal.index:
            m = massbal.loc[wy]
            # the mass-balance CSV's Season_Complete flag is the
            # admission decision (it includes SEASON_OVERRIDE_WYS);
            # fall back to the numeric screen for older CSVs
            if "Season_Complete" in m.index:
                season_ok = bool(m["Season_Complete"])
            else:
                season_ok = m.get("flood_season_missing", np.inf) \
                    <= MAX_SEASON_MISSING_DAYS
            if season_ok:
                for c in ("Three_Day", "Five_Day"):
                    if np.isfinite(m.get(c, np.nan)):
                        r[c] = m[c]
                        r["Durations_Source"] = "daily_massbalance"
                # One_day fallback: only when the hourly record didn't
                # supply it -- one-day max of the daily unreg
                if not np.isfinite(r["One_day"]) \
                        and np.isfinite(m.get("One_day", np.nan)):
                    r["One_day"] = m["One_day"]
                    r["One_day_Source"] = "daily_massbalance"
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

    print(f"Reading pre-reg peaks {PRE_REG_PEAKS_CSV}")
    pre_pk = prereg_peaks()
    print(f"  {len(pre_pk)} pre-reg USGS peaks (WY<= {PRE_REG_LAST_WY})")
    print(f"Reading daily record for pre-reg durations from {DAILY_DSS}")
    daily = read_daily()
    pre_dur = prereg_durations(daily)
    print(f"  {len(pre_dur)} pre-reg WYs with complete-season durations")

    table = assemble_record(peaks, estimates, massbal, pre_pk, pre_dur)
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
