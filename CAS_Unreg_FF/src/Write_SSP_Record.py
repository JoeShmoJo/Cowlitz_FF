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
diag_dir = os.path.join(PROJECT_DIR, "output", "diagnostics")
OUT_FLAGS_CSV = os.path.join(diag_dir, "record_qa_flags.csv")
OUT_ADJ_CSV = os.path.join(diag_dir, "record_monotonic_adjustments.csv")
OUT_EXCL_CSV = os.path.join(diag_dir, "record_excluded_wys.csv")

# Enforce Peak >= 1-day >= 3-day >= 5-day by working BACKWARDS from the
# longest duration: Five_Day is the anchor and is never changed; each
# shorter duration is raised to the longer one where the longer is
# higher, and the raise cascades upward (a lifted Three_Day can in turn
# lift One_day, which can lift Peak). Adjustments are logged to
# OUT_ADJ_CSV and summarized per WY in the record's
# Monotonic_Adjustment column, with pre-adjustment values preserved in
# the *_Raw columns. Set False to write the record unadjusted (the QA
# flags file is produced either way).
ENFORCE_MONOTONIC = True

# Water years to drop from the record ENTIRELY, whatever source would
# have supplied them (hourly holdout, regression fill, mass balance,
# pre-reg USGS). Use this to reject a year on judgment -- e.g. a year
# whose hourly record covers the regulated peak event (so the event
# screen passes) but whose season coverage is so low that a larger
# unregulated event could be hiding in the missing part.
# SEASON_OVERRIDE_WYS in Unreg_Durations_MassBalance.py does NOT do
# this: it governs only the daily mass-balance durations, so removing a
# WY there still leaves its Peak and hourly One_day in the record.
# Give a reason for each -- it is written to the exclusions log.
EXCLUDE_WYS = {
    # 2001: "39% Oct-Mar unreg coverage; annual max may fall in the gap",
}

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

# Screen for the hourly-derived peaks.
# Season-wide coverage is the WRONG test: WY2016 misses 25% of Oct-Mar
# yet covers the annual flood completely (nearest gap 43 days away), so
# a 0.9 season screen threw away a directly computed 135,790 cfs peak
# in favor of a regression estimate.
# Peak-timing offset is ALSO the wrong test for record acceptance: the
# unregulated annual maximum may legitimately occur on a different
# storm than the regulated annual maximum, because regulation decides
# which event yields the biggest regulated flow (WY2018 and WY2021 have
# 100% coverage yet offsets of 857 and -249 hours).
# The correct test: does the unregulated record actually have data
# covering the basin's biggest event, as dated by the regulated peak?
# That is unreg_cov_at_reg_peak from WY_Peak_Records.py.
MIN_EVENT_COVERAGE = 0.9
# Fallback for CSVs written before that column existed:
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
    excluded = []
    for wy in years:
        r = {"WY": wy, "Peak": np.nan, "Peak_Source": "",
             "One_day": np.nan, "One_day_Source": "",
             "Three_Day": np.nan, "Five_Day": np.nan,
             "Durations_Source": "",
             # --- audit / paper trail: both candidates and the reason ---
             "Peak_Screen": "", "Peak_Hourly_Candidate": np.nan,
             "Peak_Regression_Candidate": np.nan,
             "Peak_Offset_hrs": np.nan, "Unreg_Coverage_OctMar": np.nan,
             "Unreg_Cov_At_Reg_Peak": np.nan}

        # --- Peak and One_day from the hourly record (same-storm test) ---
        if wy in peaks.index:
            p = peaks.loc[wy]
            hourly_pk = p.get("unreg_peak_1hr", np.nan)
            evcov = p.get("unreg_cov_at_reg_peak", np.nan)
            seascov = p.get("unreg_coverage_octmar", np.nan)
            r["Peak_Hourly_Candidate"] = hourly_pk
            r["Peak_Offset_hrs"] = p.get("peak_offset_hrs", np.nan)
            r["Unreg_Coverage_OctMar"] = seascov
            r["Unreg_Cov_At_Reg_Peak"] = evcov
            if np.isfinite(evcov):
                event_ok = evcov >= MIN_EVENT_COVERAGE
                basis = f"event coverage {evcov:.2f}"
            else:   # legacy CSV without the event-coverage column
                event_ok = np.isfinite(seascov) \
                    and seascov >= MIN_UNREG_COVERAGE
                basis = f"season coverage {seascov:.2f} (legacy screen)"
            if np.isfinite(hourly_pk) and event_ok:
                r["Peak"] = hourly_pk
                r["Peak_Source"] = "hourly_holdout"
                r["Peak_Screen"] = (
                    f"hourly accepted: unreg record covers the regulated "
                    f"peak event ({basis})")
                if np.isfinite(p.get("unreg_peak_1day", np.nan)):
                    r["One_day"] = p["unreg_peak_1day"]
                    r["One_day_Source"] = "hourly_holdout"
            elif np.isfinite(hourly_pk):
                r["Peak_Screen"] = (
                    f"hourly REJECTED: unreg record does not cover the "
                    f"regulated peak event ({basis})")
            else:
                r["Peak_Screen"] = "no hourly unreg peak computed"

        # --- Peak gap fill from the adopted regression ---
        if wy in estimates.index:
            e = estimates.loc[wy]
            r["Peak_Regression_Candidate"] = e.get("unreg_peak_est", np.nan)
        if not np.isfinite(r["Peak"]) and wy in estimates.index:
            if np.isfinite(r["Peak_Regression_Candidate"]):
                r["Peak"] = r["Peak_Regression_Candidate"]
                r["Peak_Source"] = "dS2day_regression"
                if not r["Peak_Screen"]:
                    r["Peak_Screen"] = "no hourly record -> regression"
                else:
                    r["Peak_Screen"] += " -> regression"

        # --- pre-regulation WYs: observed USGS peaks + daily durations ---
        if wy <= PRE_REG_LAST_WY:
            if wy in prereg_pk.index and np.isfinite(prereg_pk[wy]):
                r["Peak"] = prereg_pk[wy]
                r["Peak_Source"] = "usgs_peak_prereg"
                r["Peak_Screen"] = "pre-regulation: USGS peak record"
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
    df = df[df[["Peak", "One_day", "Three_Day", "Five_Day"]]
            .notna().any(axis=1)]
    # --- explicit whole-WY exclusions (logged, not silently dropped) ---
    drop = [wy for wy in EXCLUDE_WYS if wy in df.index]
    if drop:
        ex = df.loc[drop, ["Peak", "Peak_Source", "One_day",
                           "One_day_Source", "Three_Day", "Five_Day",
                           "Durations_Source"]].copy()
        ex.insert(0, "exclusion_reason",
                  [EXCLUDE_WYS[wy] for wy in drop])
        os.makedirs(os.path.dirname(OUT_EXCL_CSV), exist_ok=True)
        ex.to_csv(OUT_EXCL_CSV)
        print(f"\nEXCLUDED {len(drop)} WY(s) entirely -> {OUT_EXCL_CSV}")
        print(ex[["exclusion_reason", "Peak", "Peak_Source",
                  "One_day"]].to_string())
        df = df.drop(index=drop)
    return df


def check_monotonicity(table):
    """Flag physically impossible duration ordering.

    A longer-duration average can never exceed a shorter-duration one
    drawn from the SAME series: the shorter window is free to sit on
    the wettest part of the longer one. Violations therefore mean one
    of two things, and the flag says which:

      cross_source  the two values come from different series or
                    methods (e.g. Peak from the USGS instantaneous
                    record vs One_day from the hourly rolling mean, or
                    One_day from hourly vs Three_Day from the daily
                    mass balance). Explainable, but worth reviewing --
                    it means the two durations disagree about the event.
      SAME_SOURCE   both values came from the same series. This should
                    be impossible and indicates a computation or data
                    problem; investigate before using the record.

    Returns a DataFrame of violations (empty if none).
    """
    src_col = {"Peak": "Peak_Source", "One_day": "One_day_Source",
               "Three_Day": "Durations_Source",
               "Five_Day": "Durations_Source"}
    pairs = [("Peak", "One_day"), ("One_day", "Three_Day"),
             ("Three_Day", "Five_Day"),
             # non-adjacent pairs catch cases the adjacent ones miss
             ("Peak", "Three_Day"), ("Peak", "Five_Day"),
             ("One_day", "Five_Day")]
    rows = []
    for wy, r in table.iterrows():
        for short, long_ in pairs:
            a, b = r.get(short, np.nan), r.get(long_, np.nan)
            if not (np.isfinite(a) and np.isfinite(b)) or b <= a:
                continue
            sa = str(r.get(src_col[short], "") or "")
            sb = str(r.get(src_col[long_], "") or "")
            rows.append({
                "WY": wy,
                "violation": f"{long_} > {short}",
                "shorter_duration": short, "shorter_value": a,
                "shorter_source": sa,
                "longer_duration": long_, "longer_value": b,
                "longer_source": sb,
                "excess_cfs": b - a,
                "excess_pct": 100.0 * (b - a) / a if a else np.nan,
                "flag": "SAME_SOURCE" if sa == sb else "cross_source",
            })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["flag", "WY", "violation"],
                            ascending=[True, True, True])
    return df


def enforce_monotonicity(table):
    """Raise shorter durations to match longer ones, longest first.

    Works backwards 5-day -> 3-day -> 1-day -> Peak. Five_Day is the
    anchor and is never modified. Each step compares a duration with the
    next-longer one AFTER that longer one has been finalized, so a raise
    cascades upward through the shorter durations.

    Only values that already exist are adjusted -- a missing shorter
    duration is left missing rather than fabricated from a longer one.

    Returns (adjusted_table, adjustments_log).
    """
    t = table.copy()
    order = ["Five_Day", "Three_Day", "One_day", "Peak"]  # longest first
    for col in order[1:]:
        t[f"{col}_Raw"] = t[col]
    t["Monotonic_Adjustment"] = ""
    src_col = {"Peak": "Peak_Source", "One_day": "One_day_Source",
               "Three_Day": "Durations_Source"}
    log = []
    for wy, r in t.iterrows():
        notes = []
        for longer, shorter in zip(order[:-1], order[1:]):
            lv = t.at[wy, longer]
            sv = t.at[wy, shorter]
            if not (np.isfinite(lv) and np.isfinite(sv)) or sv >= lv:
                continue
            t.at[wy, shorter] = lv
            log.append({
                "WY": wy, "duration": shorter,
                "original_value": sv, "adjusted_value": lv,
                "raised_to_match": longer,
                "delta_cfs": lv - sv,
                "delta_pct": 100.0 * (lv - sv) / sv if sv else np.nan,
                "original_source": str(r.get(src_col.get(shorter, ""), "")),
                "matched_source": str(r.get(src_col.get(longer, ""), "")),
            })
            notes.append(f"{shorter} raised {sv:,.0f} -> {lv:,.0f} "
                         f"to match {longer}")
        if notes:
            t.at[wy, "Monotonic_Adjustment"] = "; ".join(notes)
    return t, pd.DataFrame(log)


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
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    table.to_csv(OUT_CSV)
    print(f"\nAssembled {len(table)} WYs -> {OUT_CSV}")
    for col in ("Peak", "One_day", "Three_Day", "Five_Day"):
        n = int(table[col].notna().sum())
        print(f"  {col}: {n} WYs")
    n_reg = int((table["Peak_Source"] == "dS2day_regression").sum())
    print(f"  ({n_reg} peaks from the dS_2day regression fill)")

    # --- paper trail: peaks rejected by the same-storm screen ---
    rejected = table[table["Peak_Screen"].str.contains("REJECTED",
                                                      na=False)]
    if len(rejected):
        print(f"\n{len(rejected)} WY(s) where the hourly peak was rejected "
              "(flood event fell in a gap); regression used instead:")
        print(rejected[["Peak_Hourly_Candidate", "Peak_Offset_hrs",
                        "Peak"]].to_string())
        higher = rejected[rejected["Peak_Hourly_Candidate"]
                          > rejected["Peak"]]
        if len(higher):
            print("  NOTE: for WY "
                  f"{list(higher.index)} the REJECTED hourly value exceeds "
                  "the adopted regression estimate. The hourly value is an "
                  "observed event and thus a valid lower bound -- review "
                  "whether it should be adopted instead.")

    # --- advisory: hourly peaks accepted despite thin season coverage ---
    thin = table[(table["Peak_Source"] == "hourly_holdout")
                 & (table["Unreg_Coverage_OctMar"] < MIN_UNREG_COVERAGE)]
    if len(thin):
        print(f"\nADVISORY: {len(thin)} hourly peak(s) accepted because the "
              "unreg record covers the regulated peak event, but season "
              f"coverage is below {MIN_UNREG_COVERAGE} -- a larger "
              "unregulated event could lie in the uncovered part. Review; "
              "use EXCLUDE_WYS to reject any you don't trust:")
        print(thin[["Peak", "Unreg_Coverage_OctMar",
                    "Unreg_Cov_At_Reg_Peak", "Peak_Offset_hrs"]].to_string())

    # --- QA: duration ordering (BEFORE any adjustment) ---
    flags = check_monotonicity(table)
    os.makedirs(os.path.dirname(OUT_FLAGS_CSV), exist_ok=True)
    flags.to_csv(OUT_FLAGS_CSV, index=False)
    if len(flags):
        same = flags[flags["flag"] == "SAME_SOURCE"]
        cross = flags[flags["flag"] == "cross_source"]
        print(f"\nQA: {len(flags)} duration-ordering violation(s) "
              f"-> {OUT_FLAGS_CSV}")
        if len(same):
            print(f"  *** {len(same)} SAME_SOURCE violation(s) -- should be "
                  "impossible; investigate before using the record:")
            print(same[["WY", "violation", "shorter_value", "longer_value",
                        "excess_pct"]].to_string(index=False))
        if len(cross):
            print(f"  {len(cross)} cross_source violation(s) (durations from "
                  "different series disagree about the event):")
            print(cross[["WY", "violation", "shorter_source",
                         "longer_source", "excess_pct"]].to_string(
                             index=False))
    else:
        print("\nQA: duration ordering OK "
              "(Peak >= 1-day >= 3-day >= 5-day everywhere)")

    # --- enforce monotonicity, working backwards from 5-day ---
    if ENFORCE_MONOTONIC:
        table, adj = enforce_monotonicity(table)
        adj.to_csv(OUT_ADJ_CSV, index=False)
        if len(adj):
            print(f"\nMonotonic enforcement: {len(adj)} value(s) raised "
                  f"across {adj.WY.nunique()} WY(s) -> {OUT_ADJ_CSV}")
            print(adj[["WY", "duration", "original_value",
                       "adjusted_value", "raised_to_match",
                       "delta_pct"]].to_string(index=False))
        else:
            print("\nMonotonic enforcement: no adjustments needed.")
        residual = check_monotonicity(table)
        if len(residual):
            print("  *** WARNING: ordering violations remain after "
                  "enforcement -- inspect enforce_monotonicity():")
            print(residual[["WY", "violation"]].to_string(index=False))
        table.to_csv(OUT_CSV)   # rewrite with adjusted values + audit
        print(f"  record rewritten with adjustments -> {OUT_CSV}")

    print(f"\nWriting {OUT_DSS}")
    write_dss(table)
    print("\nDone. In HEC-SSP, import the *-MOS-HOLDOUT pathnames.")


if __name__ == "__main__":
    main()
