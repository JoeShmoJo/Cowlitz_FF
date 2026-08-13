#Adjusted_Peak_Record.py
# -*- coding: utf-8 -*-
"""
Build an adjusted Castle Rock annual peak record from the two ResSim runs.

PURPOSE
-------
Two ResSim runs share the same observed hydrology and the same current
operating rules, and differ only in the reservoir's starting pool:

    WCM_RC   pool starts at the WCM rule curve      (what the rules assume)
    Obs_RC   pool starts at the OBSERVED elevation  (what was actually there)

Their difference at Castle Rock therefore isolates the value of the starting
pool alone. Where WCM_RC > Obs_RC, starting at the rule curve produced a HIGHER
regulated peak than starting from the observed pool -- that is, the observed
pool held more storage than the rule curve assumes, and the historical record
benefited from storage the rules do not credit.

The adjustment removes that benefit from the observed record:

    adjusted_peak = usgs_peak + (wcm_peak - obs_peak)     when wcm > obs
    adjusted_peak = usgs_peak                             otherwise

The result is an observed-magnitude peak record placed on a consistent
rule-curve starting-pool basis, suitable for frequency analysis.

WHY THE SCREENING MATTERS
-------------------------
The difference is only meaningful if all three peaks describe the SAME storm.
If the USGS peak came from a May snowmelt event while the simulated peaks came
from a December rain event, their difference is not an operational effect and
adding it to the observed peak is meaningless.

Two independent screens are applied, and a year must pass both:

  1. TIMING. All three peak times must fall within EVENT_WINDOW_DAYS of each
     other. This is the primary screen.
  2. CONTAINMENT. The USGS peak date must fall inside the Obs_RC simulation
     window for that year (from the mapping CSV). The Obs_RC run only simulates
     a 31-day window around each year's WCM_RC peak, so if the observed peak
     falls outside that window, the Obs_RC "annual peak" is not the same event
     by construction -- no matter how the dates happen to line up.

THE THIRD SCREEN -- REGULATED ABOVE UNREGULATED AT A LARGE EVENT
----------------------------------------------------------------
A third screen is applied after the adjustment, because it can only be
evaluated once the adjusted peak exists:

  3. PHYSICALITY. If the adjusted regulated peak sits ABOVE the unregulated
     peak for the same event AND is at or above REG_OVER_UNREG_THRESHOLD_CFS,
     the year is screened out. A reservoir cannot make a flood bigger, so at a
     large event this combination means something in the chain is wrong -- the
     unregulated record, the ResSim runs, or the adjustment itself -- and the
     year cannot be used to describe the unreg-reg relationship. Below the
     threshold the crossing is expected and is NOT screened: minimum releases
     and refill drawdown legitimately put more water in the river than nature
     would, and the record is wanted for large events anyway.

     The threshold is a user setting. It typically screens only a couple of
     years, but those are exactly the years that would distort the upper end of
     the unreg-reg relationship, which is the part being relied on.

SCREENED YEARS DO NOT GO DOWNSTREAM
-----------------------------------
`screen_passed` is the AND of all three screens, and it is what
`#Critical_Duration_Adjusted.py` filters on, so a screened year reaches neither
the critical duration fits nor the unreg-reg scatter and frequency plots.
Screened years stay in `adjusted_peaks.csv` with their values, their
`screen_code` and the reason -- nothing is silently dropped -- and are listed
again on their own in `adjusted_peaks_screened_out.csv`. With
OMIT_SCREENED_FROM_EXPORTS they are held out of the SSP CSV and the DSS record,
which are the products a later step consumes. Set DROP_FAILED_SCREEN = True to
drop them from `adjusted_peaks.csv` as well (not recommended -- it destroys the
paper trail).

KNOWN CASE -- WY1980
--------------------
The observed 97,000 cfs peak on 18 May 1980 is the Mount St. Helens eruption
(18 May 1980), a debris-flow event with no meteorological analogue. Both
simulations peak on 18 Dec 1979 at 38,587 cfs from an ordinary winter storm.
This year fails the timing screen by 152 days and is never adjusted. It is
flagged explicitly in the output so the reason is recorded rather than inferred.

Whether WY1980 belongs in a flood frequency analysis at all is a separate
question this script does not decide.

REGULATED vs UNREGULATED CHECK
------------------------------
The reservoir cannot make a flood bigger, so the adjusted regulated peak must
sit below the unregulated peak for the same event. A regulated peak above the
unregulated peak is only defensible at the extreme low end, where a minimum
release or a refill drawdown puts more water in the river than nature would.
Every year is therefore compared against the unregulated peak from the ResSim
unregulated period-of-record run -- both the peak of the SAME event (within
UNREG_WINDOW_DAYS of the observed peak) and the water-year unregulated maximum,
so a timing mismatch cannot hide a violation. The check is reported and recorded
in the CSV; it does not change any value. If the unregulated run is not
available, UNREG_FALLBACK_CSV (the adopted unregulated record from
CAS_Unreg_FF) is used instead and the source is recorded per year.

OUTPUTS
-------
  adjusted_peaks.csv        every shared year, all three peaks, the adjustment,
                            screening results, the reason for each decision, and
                            the unregulated comparison
  adjusted_peaks_screened_out.csv
                            the omitted years only, with the screen that caught
                            each one -- the documentation of what was left out
  adjusted_peaks_ssp.csv    WY and adjusted peak only, for HEC-SSP import
                            (screened years omitted)
  adjusted_peaks.png        comparison and adjustment magnitude plots
  event_screening.png       timing spread per year, with the screen threshold
  //CASTLE ROCK/FLOW-ANNUAL PEAK-ADJUSTED//IR-CENTURY/<F>/  written to DSS
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
WCM_DSS = r"../output/ResSim_WCM_RC.dss"
OBS_DSS = r"../output/ResSim_Obs_RC.dss"
RESULT_DSS_VERSION = 6
PATH_WCM = "//CASTLEROCK_NWS/FLOW/*/1Hour/ResSim_WCM_RC/"
PATH_OBS = "//CASTLEROCK_NWS/FLOW/*/1Hour/ResSim_Obs_RC/"

USGS_PEAKS_CSV = r"../../CAS_Unreg_FF/data/CastleRock_USGS_peaks.csv"
OBS_MAPPING_CSV = r"../output/ensemble_obs_rc_mapping.csv"

OUT_CSV = r"../output/adjusted_peaks.csv"
OUT_SCREENED_CSV = r"../output/adjusted_peaks_screened_out.csv"
OUT_SSP_CSV = r"../output/adjusted_peaks_ssp.csv"
OUT_DSS = r"../output/adjusted_peaks.dss"
OUT_DSS_VERSION = 6
# NOTE: written as a regular 1DAY record, NOT IR-CENTURY. Handing an IR-CENTURY
# pathname to put_ts with a regular container SEGFAULTS this pydsstools build --
# it does not raise. The peak sits on its observed peak date; every other day is
# the missing sentinel, so it plots as points in DSSVue.
OUT_DSS_PATH = "/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW-ANNUAL PEAK-ADJUSTED//1DAY/ADJUSTED/"
PLOT_STEM = r"../output/diagnostics/adjusted_peaks"

WATER_YEAR_START_MONTH = 10

# --- how the simulated peaks are picked ---
#   "event"  : the simulated peak WITHIN +/- EVENT_WINDOW_DAYS of the observed
#              peak date. This is the right choice: the two runs only differ in
#              starting pool, so the comparison must be made on the same storm.
#              The Obs_RC run simulates a 31-day window, and a low starting pool
#              can attenuate the first storm enough that a LATER storm in the
#              same window becomes that run's annual maximum. Comparing annual
#              maxima would then difference two different storms.
#   "annual" : each run's water-year maximum, whenever it occurs. Kept for
#              comparison; it screens out roughly twice as many years.
PEAK_METHOD = "event"

# --- screening ---
# All three peak times must fall within this many days of each other. Also the
# half-width of the event window when PEAK_METHOD = "event".
EVENT_WINDOW_DAYS = 3
# The USGS peak must also fall inside the Obs_RC simulation window for the year.
REQUIRE_IN_OBS_WINDOW = True
# Years failing a screen: False = keep the row in adjusted_peaks.csv flagged and
# unadjusted (recommended -- it is the paper trail), True = drop it from that
# file too. Either way a screened year is excluded from everything downstream,
# because #Critical_Duration_Adjusted.py filters on screen_passed.
DROP_FAILED_SCREEN = False
# Hold screened years out of the products a later step reads: the SSP CSV and
# the DSS record. adjusted_peaks.csv always keeps them (unless DROP_FAILED_SCREEN).
OMIT_SCREENED_FROM_EXPORTS = True

# --- screen 3: regulated peak above unregulated at a large event ------------
# Applied AFTER the adjustment, since it is the ADJUSTED peak that has to sit
# below the unregulated peak. A year is screened out when BOTH hold:
#     adjusted peak  >  unregulated peak for the same event (by UNREG_TOL_CFS)
#     adjusted peak  >= REG_OVER_UNREG_THRESHOLD_CFS
# Below the threshold the crossing is expected -- minimum release and refill
# drawdown put more water in the river than nature would -- so those years are
# reported but kept. This is the setting to move if the screen is catching too
# much or too little; it is deliberately separate from UNREG_LOW_FLOW_CFS,
# which only controls how the crossing is REPORTED.
SCREEN_REG_OVER_UNREG = True
REG_OVER_UNREG_THRESHOLD_CFS = 60000.0
# Test the threshold against the adjusted regulated peak ("reg") or against the
# unregulated peak it is being compared with ("unreg"). "reg" matches the way
# the screen is described: screen where the REGULATED peak is large.
REG_OVER_UNREG_THRESHOLD_ON = "reg"
# Years whose USGS peak is below this are reported separately -- low-flow years
# are where spurious timing mismatches cluster, because the annual maximum is
# not a distinct storm.
LOW_PEAK_CFS = 20000.0
# Differences smaller than this are not worth recording as an adjustment.
MIN_ADJUSTMENT_CFS = 1.0

# --- regulated vs unregulated check -----------------------------------------
CHECK_UNREG = True
# EXTERNAL: the ResSim unregulated period-of-record run (not in the repository).
# "Ensemble--0" in the F-part is a ResSim artifact; this is a single POR run.
# Same source #Critical_Duration_Adjusted.py uses, so the two agree by
# construction.
UNREG_DSS = (r"C:\Projects\2026_Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis"
             r"\watershed\NWP_CowlitzLewis_ResSim4\rss\Unreg_POR_FIS\simulation.dss")
UNREG_PATH = "//CastleRock_NWS/Flow-UNREG/*/1Hour/Ensemble--0/"
# Used when UNREG_DSS is not reachable: the adopted unregulated record from the
# flow frequency study. Annual values only -- no event matching is possible, so
# only the water-year comparison is made and the source is recorded per year.
UNREG_FALLBACK_CSV = r"../../CAS_Unreg_FF/output/wy_record_ssp.csv"
UNREG_FALLBACK_COLS = ("WY", "Peak")
# The unregulated peak for the SAME event: the unregulated maximum within this
# many days of the observed peak. Wider than EVENT_WINDOW_DAYS on purpose -- the
# unregulated peak leads the regulated one, and a too-tight window would report
# a false violation instead of a real one.
UNREG_WINDOW_DAYS = 5
# Ignore crossings below this; DSS is single precision and the two records take
# different routes through ResSim.
UNREG_TOL_CFS = 1.0
# Above this the crossing cannot be explained by minimum releases or refill
# drawdown, and is called out as a hard failure rather than a note.
UNREG_LOW_FLOW_CFS = 20000.0

# Years never to adjust regardless of screening, with the reason recorded.
# 1980: Mount St. Helens eruption, 18 May 1980 -- not a meteorological event.
MANUAL_EXCLUSIONS = {
    1980: "Mount St. Helens eruption 18 May 1980; not a meteorological event",
}

SENTINEL = -901.0

# ----------------------------------------------------------------------------


def first_stamp(ts):
    """First timestamp of a DSS series, across pydsstools versions."""
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    # DSS uses midnight-as-2400: "01Oct1973 24:00:00" means 02Oct1973 00:00
    roll_day = False
    if " 24:" in text or text.endswith(" 2400"):
        text = text.replace(" 24:", " 00:").replace(" 2400", " 0000")
        roll_day = True
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M", "%d%b%Y %H%M%S", "%d%b%Y %H%M",
                "%d %B %Y %H:%M:%S", "%d %B %Y %H:%M"):
        try:
            stamp = pd.Timestamp(datetime.strptime(text, fmt))
            return stamp + pd.Timedelta(days=1) if roll_day else stamp
        except ValueError:
            continue
    stamp = pd.Timestamp(text)
    return stamp + pd.Timedelta(days=1) if roll_day else stamp


def series_step(ts, pathname):
    """Time step of a DSS regular series. ts.interval is in seconds."""
    seconds = int(getattr(ts, "interval", 0) or 0)
    if seconds > 0:
        return pd.Timedelta(seconds=seconds)
    e_part = pathname.split("/")[5].upper()
    lookup = {"1MIN": "1min", "15MIN": "15min", "30MIN": "30min",
              "1HOUR": "1h", "6HOUR": "6h", "12HOUR": "12h", "1DAY": "1D"}
    return pd.Timedelta(lookup.get(e_part, "1h"))


def dss_version(path):
    """DSS file version from the header: byte 12 is 6 for v6, 0 for v7.

    The ResSim runs are a mix -- the newer simulation.dss files are v7 while the
    reassembled result files here are v6 -- and pydsstools needs the version
    passed explicitly on Linux, so detect rather than assume.
    """
    with open(path, "rb") as handle:
        head = handle.read(16)
    if len(head) < 16 or head[:4] != b"ZDSS":
        return None
    return 6 if head[12] == 6 else 7


def read_dss_series(dss_file, pathname, version):
    """Read a DSS regular time series into a Series on period-BEGINNING labels."""
    dss = HecDss.Open(dss_file, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        nodata = np.array(ts.nodata, dtype=bool)
        values[nodata] = np.nan
        values[values <= -900.0] = np.nan
        step = series_step(ts, pathname)
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index().dropna()


def water_year(stamp):
    return stamp.year + (1 if stamp.month >= WATER_YEAR_START_MONTH else 0)


def annual_peaks(series):
    """Water-year maximum and its timestamp."""
    keys = [water_year(t) for t in series.index]
    grouped = series.groupby(keys)
    return pd.DataFrame({"peak": grouped.max(), "time": grouped.idxmax()})


def event_peak(series, when, half_width_days):
    """Peak within +/- half_width_days of a given time, or (nan, NaT) if absent.

    Returns the simulated response to the SAME storm the observed peak came
    from, which is what the WCM/Obs difference is supposed to isolate.
    """
    a = when - pd.Timedelta(days=half_width_days)
    b = when + pd.Timedelta(days=half_width_days)
    window = series.loc[a:b]
    if len(window) == 0:
        return np.nan, pd.NaT
    return float(window.max()), pd.Timestamp(window.idxmax())


def load_usgs_peaks(csv_path):
    """USGS annual peaks: columns Peak_Date, WY, Peak_cfs."""
    table = pd.read_csv(csv_path, parse_dates=["Peak_Date"])
    table.columns = [c.strip() for c in table.columns]
    return table.set_index("WY")


def load_obs_windows(csv_path):
    """Obs_RC simulation window per water year, for the containment screen."""
    if not os.path.isfile(csv_path):
        return None
    table = pd.read_csv(csv_path, parse_dates=["real_start", "real_end"])
    return {int(r["water_year"]): (r["real_start"], r["real_end"])
            for _, r in table.iterrows()}


def build_container(pathname, values, start_time, units, data_type, interval_min):
    """TimeSeriesContainer built the way the installed pydsstools accepts."""
    try:
        tsc = TimeSeriesContainer(pathname, len(values), interval_min,
                                  values=list(values), start_time=start_time,
                                  data_units=units, data_type=data_type)
    except TypeError:
        tsc = TimeSeriesContainer()
        tsc.pathname = pathname
        tsc.startDateTime = start_time
        tsc.numberValues = len(values)
        tsc.interval = interval_min
        tsc.values = list(values)
        tsc.units = units
        tsc.type = data_type
    return tsc


def load_unreg():
    """Unregulated Castle Rock flow for the check. Returns (kind, data).

    kind is "hourly" (a ResSim POR series, event matching possible), "annual"
    (the adopted CAS_Unreg_FF record, water-year comparison only) or None.
    """
    if not CHECK_UNREG:
        return None, None
    if os.path.isfile(UNREG_DSS):
        try:
            series = read_dss_series(UNREG_DSS, UNREG_PATH, dss_version(UNREG_DSS))
            print("Unreg   ResSim unregulated POR run, %s .. %s   %d hours"
                  % (series.index[0].date(), series.index[-1].date(), len(series)))
            return "hourly", series
        except Exception as exc:
            print("Unreg   could not read %s -- %s" % (UNREG_DSS, exc))
    else:
        print("Unreg   not found: %s" % UNREG_DSS)
    if os.path.isfile(UNREG_FALLBACK_CSV):
        wy_col, peak_col = UNREG_FALLBACK_COLS
        table = pd.read_csv(UNREG_FALLBACK_CSV)
        data = {int(row[wy_col]): float(row[peak_col])
                for _, row in table.iterrows() if pd.notna(row[peak_col])}
        print("Unreg   FALLBACK to the adopted annual record %s (%d years)"
              % (UNREG_FALLBACK_CSV, len(data)))
        print("        annual values only -- the same-event comparison is skipped")
        return "annual", data
    print("Unreg   no unregulated source available; the check is skipped")
    return None, None


def unreg_for_year(kind, data, wy_peaks, wy, t_usgs):
    """Unregulated peak for one year: same-event value and water-year maximum."""
    if kind == "hourly":
        p_evt, t_evt = event_peak(data, t_usgs, UNREG_WINDOW_DAYS)
        if wy_peaks is not None and wy in wy_peaks.index:
            p_wy = float(wy_peaks.loc[wy, "peak"])
            t_wy = pd.Timestamp(wy_peaks.loc[wy, "time"])
        else:
            p_wy, t_wy = np.nan, pd.NaT
        return p_evt, t_evt, p_wy, t_wy, "ressim_unreg_por"
    if kind == "annual":
        return np.nan, pd.NaT, float(data.get(wy, np.nan)), pd.NaT, "cas_unreg_record"
    return np.nan, pd.NaT, np.nan, pd.NaT, "none"


def unreg_columns(p_adj, p_usgs, p_evt, t_evt, p_wy, t_wy, source):
    """Compare a regulated peak against the unregulated peak. Reports only.

    The same-event unregulated peak is the reference where it exists, because it
    is the like-for-like comparison. The water-year maximum is carried alongside
    it as the timing-independent backstop: if the event value is missing, or if
    the adjusted peak clears the event value but not the year's unregulated
    maximum, that shows up here rather than being lost.
    """
    ref = p_evt if np.isfinite(p_evt) else p_wy
    over_adj = p_adj - ref if np.isfinite(ref) else np.nan
    over_usgs = p_usgs - ref if np.isfinite(ref) else np.nan
    flag = bool(np.isfinite(over_adj) and over_adj > UNREG_TOL_CFS)
    hard = bool(flag and ref >= UNREG_LOW_FLOW_CFS)
    if not np.isfinite(ref):
        note = "no unregulated value for this year"
    elif not flag:
        note = "ok: adjusted %.0f <= unreg %.0f" % (p_adj, ref)
    elif hard:
        note = ("*** adjusted %.0f EXCEEDS unreg %.0f by %.0f cfs (%.1f%%) at a "
                "flood peak -- the reservoir cannot raise a flood, so either the "
                "adjustment, the ResSim runs or the unregulated record is wrong"
                % (p_adj, ref, over_adj, 100.0 * over_adj / ref))
    else:
        note = ("adjusted %.0f exceeds unreg %.0f by %.0f cfs, but the unreg peak "
                "is below %.0f cfs -- consistent with minimum release or refill "
                "drawdown" % (p_adj, ref, over_adj, UNREG_LOW_FLOW_CFS))
    return {
        "unreg_event": p_evt,
        "t_unreg_event": t_evt.date() if pd.notna(t_evt) else pd.NaT,
        "unreg_wy_max": p_wy,
        "t_unreg_wy": t_wy.date() if pd.notna(t_wy) else pd.NaT,
        "unreg_ref": ref,
        "unreg_source": source,
        "adj_minus_unreg": over_adj,
        "usgs_minus_unreg": over_usgs,
        "adj_over_unreg": flag,
        "adj_over_unreg_at_flood": hard,
        "usgs_over_unreg": bool(np.isfinite(over_usgs) and over_usgs > UNREG_TOL_CFS),
        "adj_unreg_ratio": p_adj / ref if np.isfinite(ref) and ref > 0 else np.nan,
        "unreg_check": note,
    }


def screen_year(wy, t_usgs, t_wcm, t_obs, obs_windows):
    """Do all three peaks describe the same event? Returns (passed, reason)."""
    if wy in MANUAL_EXCLUSIONS:
        return False, "excluded: %s" % MANUAL_EXCLUSIONS[wy]

    spread = max(abs((t_wcm - t_usgs).days),
                 abs((t_obs - t_usgs).days),
                 abs((t_wcm - t_obs).days))
    if spread > EVENT_WINDOW_DAYS:
        return False, ("peak times differ by %d days (limit %d): "
                       "USGS %s, WCM %s, Obs %s"
                       % (spread, EVENT_WINDOW_DAYS, t_usgs.date(),
                          t_wcm.date(), t_obs.date()))

    if REQUIRE_IN_OBS_WINDOW and obs_windows and wy in obs_windows:
        start, end = obs_windows[wy]
        if not (start <= t_usgs <= end):
            return False, ("USGS peak %s outside the Obs_RC window %s..%s"
                           % (t_usgs.date(), pd.Timestamp(start).date(),
                              pd.Timestamp(end).date()))
    return True, "same event"


def apply_reg_over_unreg_screen(table):
    """Screen 3, then resolve the combined screen. Adds the audit columns.

    Runs on the assembled table because the adjusted peak has to exist first.
    Adds:
        screen_reg_le_unreg  True when the year passes screen 3
        screen_passed        AND of the same-event screens and screen 3
        screen_code          ok | different_event | reg_over_unreg | both
        screen_reason        why, in words
    """
    n = len(table)
    passes = np.ones(n, dtype=bool)
    reasons = [""] * n

    if SCREEN_REG_OVER_UNREG:
        for i in range(n):
            row = table.iloc[i]
            p_adj = float(row["adjusted_peak"])
            ref = float(row["unreg_ref"]) if pd.notna(row["unreg_ref"]) else np.nan
            if not np.isfinite(ref):
                continue
            if p_adj - ref <= UNREG_TOL_CFS:
                continue
            gauge = p_adj if REG_OVER_UNREG_THRESHOLD_ON == "reg" else ref
            if gauge < REG_OVER_UNREG_THRESHOLD_CFS:
                continue
            passes[i] = False
            reasons[i] = ("regulated peak %.0f exceeds the unregulated peak %.0f "
                          "by %.0f cfs (%.1f%%) at or above the %s screening "
                          "threshold of %s cfs -- a reservoir cannot raise a "
                          "flood, so the pair is not usable"
                          % (p_adj, ref, p_adj - ref, 100.0 * (p_adj - ref) / ref,
                             REG_OVER_UNREG_THRESHOLD_ON,
                             format(int(REG_OVER_UNREG_THRESHOLD_CFS), ",")))

    table = table.copy()
    table["screen_reg_le_unreg"] = passes
    same_event = table["screen_same_event"].values.astype(bool)
    table["screen_passed"] = same_event & passes

    codes, notes = [], []
    for i in range(n):
        if same_event[i] and passes[i]:
            codes.append("ok")
            notes.append("same event; adjusted peak below the unregulated peak")
        elif not same_event[i] and not passes[i]:
            codes.append("both")
            notes.append("%s; and %s" % (table["decision"].iloc[i], reasons[i]))
        elif not same_event[i]:
            codes.append("different_event")
            notes.append(str(table["decision"].iloc[i]))
        else:
            codes.append("reg_over_unreg")
            notes.append(reasons[i])
    table["screen_code"] = codes
    table["screen_reason"] = notes

    # Fold the new screen into the decision text so one column still reads as
    # the record of what happened to that year.
    for i in range(n):
        if not passes[i]:
            table.iloc[i, table.columns.get_loc("decision")] = (
                "SCREENED OUT: %s" % reasons[i])
    return table


def report_screening(table):
    """Print, and return, what was screened out and why."""
    screened = table[~table["screen_passed"]]
    print("\n" + "=" * 78)
    print("EVENT SCREENING -- years omitted from everything downstream")
    print("=" * 78)
    print("   eligible : %d of %d years" % (int(table["screen_passed"].sum()),
                                            len(table)))
    if screened.empty:
        print("   nothing screened out")
        return screened
    by_code = screened["screen_code"].value_counts()
    print("   screened : %d  (%s)"
          % (len(screened), ", ".join("%s %d" % (k, v) for k, v in by_code.items())))
    print("   screen 3 threshold: %s cfs on the %s peak"
          % (format(int(REG_OVER_UNREG_THRESHOLD_CFS), ","),
             REG_OVER_UNREG_THRESHOLD_ON))
    for _, row in screened.sort_values("WY").iterrows():
        print("      WY%d  [%s]  usgs %.0f  adjusted %.0f  unreg %s"
              % (row["WY"], row["screen_code"], row["usgs"], row["adjusted_peak"],
                 format(int(row["unreg_ref"]), ",") if pd.notna(row["unreg_ref"])
                 else "n/a"))
        print("            %s" % row["screen_reason"])
    return screened


def report_unreg(table):
    """Print the regulated vs unregulated check."""
    if "unreg_ref" not in table.columns:
        return
    have = table[table["unreg_ref"].notna()]
    print("\n" + "=" * 78)
    print("REGULATED vs UNREGULATED CHECK")
    print("=" * 78)
    if have.empty:
        print("   no unregulated values available -- nothing compared")
        return
    sources = ", ".join("%s (%d)" % (k, v) for k, v in
                        have["unreg_source"].value_counts().items())
    print("   years compared : %d of %d   source: %s"
          % (len(have), len(table), sources))

    hard = have[have["adj_over_unreg_at_flood"]]
    soft = have[have["adj_over_unreg"] & ~have["adj_over_unreg_at_flood"]]
    print("   adjusted peak ABOVE unregulated at a flood peak (>= %.0f cfs): %d"
          % (UNREG_LOW_FLOW_CFS, len(hard)))
    for _, row in hard.sort_values("adj_minus_unreg", ascending=False).iterrows():
        print("      WY%d  adjusted %.0f vs unreg %.0f  (+%.0f cfs, %.2fx)"
              % (row["WY"], row["adjusted_peak"], row["unreg_ref"],
                 row["adj_minus_unreg"], row["adj_unreg_ratio"]))
        if row["adjusted"]:
            print("            the adjustment (+%.0f) is what pushes it over: "
                  "USGS %.0f vs unreg %.0f"
                  % (row["delta_wcm_minus_obs"], row["usgs"], row["unreg_ref"]))
        if bool(row["usgs_over_unreg"]):
            print("            the OBSERVED peak already exceeds the "
                  "unregulated peak -- the problem is upstream of this script")
    if len(hard):
        print("      A reservoir cannot raise a flood peak. Look at, in order:")
        print("        1. the unregulated record -- is it the same event, and "
              "does the POR run cover this year?")
        print("        2. the ResSim runs -- WCM_RC releasing more than the "
              "operation set intends")
        print("        3. the adjustment itself -- WCM_RC and Obs_RC peaking on "
              "different storms inside the window")
    print("   adjusted peak above unregulated at a low peak (< %.0f cfs): %d"
          % (UNREG_LOW_FLOW_CFS, len(soft)))
    for _, row in soft.sort_values("WY").iterrows():
        print("      WY%d  adjusted %.0f vs unreg %.0f  (+%.0f cfs)"
              % (row["WY"], row["adjusted_peak"], row["unreg_ref"],
                 row["adj_minus_unreg"]))
    if not len(hard) and not len(soft):
        print("   OK -- every adjusted peak sits at or below its unregulated peak")

    # Event vs water-year backstop: a peak that clears the same-event unreg value
    # but not the year's unreg maximum means the two records peak on different
    # storms, which is a timing problem rather than a magnitude one.
    both = have[have["unreg_event"].notna() & have["unreg_wy_max"].notna()]
    if len(both):
        mismatch = both[(both["adjusted_peak"] <= both["unreg_event"])
                        & (both["adjusted_peak"] > both["unreg_wy_max"])]
        if len(mismatch):
            print("   %d year(s) clear the same-event unreg peak but exceed the "
                  "water-year unreg maximum: %s"
                  % (len(mismatch), ", ".join("WY%d" % w for w in mismatch["WY"])))


def plot_peaks(table, stem):
    """Three peak records and the adjustment applied."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    x = table["WY"].values
    adjusted = table["adjusted"].values
    passed = table["screen_passed"].values

    ax = axes[0]
    ax.plot(x, table["usgs"], color="0.25", lw=1.4, marker="o", ms=3.5,
            label="USGS observed")
    ax.plot(x, table["wcm"], color="#c0392b", lw=1.2, marker="s", ms=3,
            label="ResSim WCM_RC (rule curve start)")
    ax.plot(x, table["obs"], color="#2c7fb8", lw=1.2, marker="^", ms=3,
            label="ResSim Obs_RC (observed pool start)")
    ax.plot(x[adjusted], table.loc[adjusted, "adjusted_peak"], color="#16a085",
            lw=0, marker="*", ms=11, label="Adjusted peak")
    if "unreg_ref" in table.columns and table["unreg_ref"].notna().any():
        ax.plot(x, table["unreg_ref"], color="#8e44ad", lw=1.0, ls="--",
                marker="v", ms=3, alpha=0.8, label="Unregulated peak (ceiling)")
        over = table["adj_over_unreg"].fillna(False).values.astype(bool)
        if over.any():
            ax.plot(x[over], table.loc[over, "adjusted_peak"], lw=0, marker="o",
                    ms=14, mfc="none", mec="#c0392b", mew=1.6,
                    label="Adjusted peak ABOVE unregulated")
    for _, row in table[~table["screen_passed"]].iterrows():
        ax.axvspan(row["WY"] - 0.4, row["WY"] + 0.4, color="0.85", alpha=0.5, zorder=0)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("Castle Rock annual peaks -- observed, both ResSim runs, and the "
                 "adjusted record\n(grey bands = year failed event screening, "
                 "carried through unadjusted)", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    delta = (table["wcm"] - table["obs"]).fillna(0.0)
    colors = ["#16a085" if a else "0.7" for a in adjusted]
    ax.bar(x, delta, color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("WCM_RC - Obs_RC (cfs)")
    ax.set_xlabel("Water year")
    ax.set_title("Adjustment applied: the amount the observed peak is increased "
                 "(green = applied, grey = not applied)", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig("%s.png" % stem, dpi=150)
    plt.close(fig)


def plot_screening(table, stem):
    """Timing spread per year against the screening threshold."""
    fig, ax = plt.subplots(figsize=(14, 5.5))
    x = table["WY"].values
    palette = {"ok": "#16a085", "different_event": "#c0392b",
               "reg_over_unreg": "#e67e22", "both": "#7d3c98"}
    colors = [palette.get(c, "#c0392b") for c in table["screen_code"]]
    ax.bar(x, table["spread_days"].fillna(0.0), color=colors)
    ax.axhline(EVENT_WINDOW_DAYS, color="k", lw=1.1, ls="--")
    ax.text(x.min(), EVENT_WINDOW_DAYS + 1, "screen: %d days" % EVENT_WINDOW_DAYS,
            fontsize=9)
    for _, row in table[~table["screen_passed"]].iterrows():
        ax.annotate("%d" % row["WY"], (row["WY"], row["spread_days"]
                                       if np.isfinite(row["spread_days"]) else 0.0),
                    xytext=(0, 3), textcoords="offset points", ha="center",
                    fontsize=7, rotation=90)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylabel("Max spread between the three peak times (days)")
    ax.set_xlabel("Water year")
    ax.set_title("Event screening -- do the USGS, WCM_RC and Obs_RC peaks "
                 "describe the same storm, and does the adjusted peak stay "
                 "below the unregulated peak?")
    handles = [Line2D([], [], color=palette["ok"], lw=6, label="Eligible"),
               Line2D([], [], color=palette["different_event"], lw=6,
                      label="Different event, screened out"),
               Line2D([], [], color=palette["reg_over_unreg"], lw=6,
                      label="Reg above unreg at >= %s cfs, screened out"
                            % format(int(REG_OVER_UNREG_THRESHOLD_CFS), ",")),
               Line2D([], [], color=palette["both"], lw=6,
                      label="Both screens failed")]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig("%s_screening.png" % stem, dpi=150)
    plt.close(fig)


def main():
    for path in (os.path.dirname(OUT_CSV), os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    wcm = read_dss_series(WCM_DSS, PATH_WCM, RESULT_DSS_VERSION)
    obs = read_dss_series(OBS_DSS, PATH_OBS, RESULT_DSS_VERSION)
    peaks_wcm = annual_peaks(wcm)
    peaks_obs = annual_peaks(obs)
    usgs = load_usgs_peaks(USGS_PEAKS_CSV)
    obs_windows = load_obs_windows(OBS_MAPPING_CSV)
    unreg_kind, unreg_data = load_unreg()
    unreg_wy_peaks = annual_peaks(unreg_data) if unreg_kind == "hourly" else None

    shared = sorted(set(peaks_wcm.index) & set(peaks_obs.index) & set(usgs.index))

    print("=" * 78)
    print("WCM_RC  %s .. %s   %d water years"
          % (wcm.index[0].date(), wcm.index[-1].date(), len(peaks_wcm)))
    print("Obs_RC  %s .. %s   %d water years"
          % (obs.index[0].date(), obs.index[-1].date(), len(peaks_obs)))
    print("USGS    %d water years (%d..%d)"
          % (len(usgs), usgs.index.min(), usgs.index.max()))
    print("Shared  %d water years (%d..%d)" % (len(shared), shared[0], shared[-1]))
    print("Method  %s peaks (%s)"
          % (PEAK_METHOD,
             "simulated peak within +/-%d days of the observed peak" % EVENT_WINDOW_DAYS
             if PEAK_METHOD == "event" else "each run's water-year maximum"))
    print("Screen  peaks within %d days; USGS peak inside the Obs_RC window: %s"
          % (EVENT_WINDOW_DAYS, REQUIRE_IN_OBS_WINDOW))
    print("=" * 78)

    rows = []
    for wy in shared:
        t_usgs = pd.Timestamp(usgs.loc[wy, "Peak_Date"])
        p_usgs = float(usgs.loc[wy, "Peak_cfs"])
        u_evt, u_t_evt, u_wy, u_t_wy, u_src = unreg_for_year(
            unreg_kind, unreg_data, unreg_wy_peaks, wy, t_usgs)
        if PEAK_METHOD == "event":
            p_wcm, t_wcm = event_peak(wcm, t_usgs, EVENT_WINDOW_DAYS)
            p_obs, t_obs = event_peak(obs, t_usgs, EVENT_WINDOW_DAYS)
            if not (np.isfinite(p_wcm) and np.isfinite(p_obs)):
                row = {
                    "WY": wy, "usgs": p_usgs, "wcm": np.nan, "obs": np.nan,
                    "delta_wcm_minus_obs": np.nan, "adjusted_peak": p_usgs,
                    "adjusted": False, "screen_same_event": False,
                    "spread_days": np.nan, "t_usgs": t_usgs.date(),
                    "t_wcm": pd.NaT, "t_obs": pd.NaT,
                    "low_peak_year": p_usgs < LOW_PEAK_CFS,
                    "decision": ("no adjustment: no simulated data within %d days "
                                 "of the observed peak %s -- the Obs_RC window "
                                 "does not cover this storm"
                                 % (EVENT_WINDOW_DAYS, t_usgs.date()))}
                row.update(unreg_columns(p_usgs, p_usgs, u_evt, u_t_evt,
                                         u_wy, u_t_wy, u_src))
                rows.append(row)
                continue
        else:
            t_wcm = pd.Timestamp(peaks_wcm.loc[wy, "time"])
            t_obs = pd.Timestamp(peaks_obs.loc[wy, "time"])
            p_wcm = float(peaks_wcm.loc[wy, "peak"])
            p_obs = float(peaks_obs.loc[wy, "peak"])
        spread = max(abs((t_wcm - t_usgs).days), abs((t_obs - t_usgs).days),
                     abs((t_wcm - t_obs).days))
        passed, reason = screen_year(wy, t_usgs, t_wcm, t_obs, obs_windows)
        delta = p_wcm - p_obs
        # ONE-SIDED BY DESIGN: the observed peak is only ever INCREASED. A
        # negative difference means the observed pool started ABOVE the rule
        # curve and the historical operation was already at least as good as
        # following the WCM; there is nothing to remove, so the peak stands.
        adjust = bool(passed and delta > MIN_ADJUSTMENT_CFS)
        if adjust:
            decision = "adjusted +%.0f cfs" % delta
        elif passed and delta < 0:
            decision = ("no adjustment: Obs_RC peak exceeds WCM_RC by %.0f cfs "
                        "-- observed pool started above the rule curve, so the "
                        "observed peak is not reduced" % abs(delta))
        elif passed:
            decision = "no adjustment: difference below %.0f cfs" % MIN_ADJUSTMENT_CFS
        else:
            decision = "no adjustment: %s" % reason
        p_adjusted = p_usgs + delta if adjust else p_usgs
        row = {
            "WY": wy,
            "usgs": p_usgs, "wcm": p_wcm, "obs": p_obs,
            "delta_wcm_minus_obs": delta,
            "adjusted_peak": p_adjusted,
            "adjusted": adjust,
            "screen_same_event": passed,
            "spread_days": spread,
            "t_usgs": t_usgs.date(), "t_wcm": t_wcm.date(), "t_obs": t_obs.date(),
            "low_peak_year": p_usgs < LOW_PEAK_CFS,
            "decision": decision,
        }
        row.update(unreg_columns(p_adjusted, p_usgs, u_evt, u_t_evt,
                                 u_wy, u_t_wy, u_src))
        rows.append(row)

    table = pd.DataFrame(rows)
    table = apply_reg_over_unreg_screen(table)
    screened = report_screening(table)
    screened[["WY", "usgs", "wcm", "obs", "delta_wcm_minus_obs", "adjusted_peak",
              "unreg_ref", "unreg_source", "adj_minus_unreg", "spread_days",
              "t_usgs", "t_wcm", "t_obs", "screen_code", "screen_reason"]].to_csv(
        OUT_SCREENED_CSV, index=False, float_format="%.1f")

    if DROP_FAILED_SCREEN:
        table = table[table["screen_passed"]].reset_index(drop=True)

    table.to_csv(OUT_CSV, index=False, float_format="%.1f")
    # Everything that feeds a later step carries eligible years only.
    export = table[table["screen_passed"]] if OMIT_SCREENED_FROM_EXPORTS else table
    export[["WY", "adjusted_peak"]].to_csv(OUT_SSP_CSV, index=False,
                                           float_format="%.0f")

    if os.path.exists(OUT_DSS):
        os.remove(OUT_DSS)
    dss = HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION)
    try:
        span = pd.date_range(pd.Timestamp(int(table["WY"].min()) - 1, 10, 1),
                             pd.Timestamp(int(table["WY"].max()), 9, 30), freq="D")
        values = pd.Series(SENTINEL, index=span)
        for _, row in export.iterrows():
            values.loc[pd.Timestamp(row["t_usgs"])] = row["adjusted_peak"]
        # DSS stamps are end-of-period: a 1DAY value is stamped at the NEXT midnight
        start_time = (span[0] + pd.Timedelta(days=1)).strftime("%d%b%Y %H:%M:%S").upper()
        dss.put_ts(build_container(OUT_DSS_PATH, values.values, start_time,
                                   "CFS", "INST-VAL", 1440))
    finally:
        dss.close()

    plot_peaks(table, PLOT_STEM)
    plot_screening(table, PLOT_STEM)

    adjusted = table[table["adjusted"]]
    failed = table[~table["screen_passed"]]
    no_gain = table[table["screen_passed"] & ~table["adjusted"]]

    print("\nADJUSTED   %d of %d years" % (len(adjusted), len(table)))
    if len(adjusted):
        print("   adjustment: median %+.0f cfs, range %+.0f to %+.0f"
              % (adjusted["delta_wcm_minus_obs"].median(),
                 adjusted["delta_wcm_minus_obs"].min(),
                 adjusted["delta_wcm_minus_obs"].max()))
        print("   as a share of the observed peak: median %+.1f%%, max %+.1f%%"
              % ((100 * adjusted["delta_wcm_minus_obs"] / adjusted["usgs"]).median(),
                 (100 * adjusted["delta_wcm_minus_obs"] / adjusted["usgs"]).max()))
    print("NOT ADJUSTED, no gain from the rule curve: %d years  "
          "(record is increase-only; peaks are never reduced)" % len(no_gain))
    print("NOT ADJUSTED, failed event screening     : %d years" % len(failed))
    for _, row in failed.iterrows():
        flag = "  [low peak %.0f cfs]" % row["usgs"] if row["low_peak_year"] else ""
        print("   WY%d  %s%s" % (row["WY"], row["decision"], flag))
    low_fail = failed[failed["low_peak_year"]]
    print("\n   of the %d screened out, %d are low-peak years (< %.0f cfs)"
          % (len(failed), len(low_fail), LOW_PEAK_CFS))

    # A cluster of adjustments at one value means the Obs_RC run is hitting an
    # operational limit rather than responding to the starting pool.
    if len(adjusted):
        rounded = adjusted["delta_wcm_minus_obs"].round(-2)
        mode = rounded.mode()
        if len(mode):
            n_mode = int((rounded == mode.iloc[0]).sum())
            if n_mode >= 3:
                print("\n   NOTE: %d of %d adjustments cluster near %+.0f cfs."
                      % (n_mode, len(adjusted), mode.iloc[0]))
                print("   A repeated identical difference usually means the two runs")
                print("   are separated by a fixed release rule, not by the starting")
                print("   pool. Worth confirming against the ResSim operation set.")
        big = adjusted[adjusted["delta_wcm_minus_obs"] > adjusted["usgs"] * 0.4]
        if len(big):
            print("\n   %d adjustments exceed 40%% of the observed peak: %s"
                  % (len(big), ", ".join("WY%d +%.0f%%" % (r["WY"],
                     100 * r["delta_wcm_minus_obs"] / r["usgs"])
                     for _, r in big.iterrows())))
    report_unreg(table)

    print("-" * 78)
    print("Peaks CSV   : %s" % OUT_CSV)
    print("Screened out: %s  (%d years, omitted downstream)"
          % (OUT_SCREENED_CSV, len(screened)))
    print("SSP CSV     : %s%s"
          % (OUT_SSP_CSV,
             "  (eligible years only)" if OMIT_SCREENED_FROM_EXPORTS else ""))
    print("DSS         : %s  %s" % (OUT_DSS, OUT_DSS_PATH))
    print("Plots       : %s.png and %s_screening.png" % (PLOT_STEM, PLOT_STEM))


main()
