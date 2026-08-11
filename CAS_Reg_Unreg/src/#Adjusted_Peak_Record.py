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

Years failing either screen are carried through UNADJUSTED, with the reason
recorded, rather than dropped. That keeps the record complete and makes every
decision auditable. Set DROP_FAILED_SCREEN = True to omit them instead.

KNOWN CASE -- WY1980
--------------------
The observed 97,000 cfs peak on 18 May 1980 is the Mount St. Helens eruption
(18 May 1980), a debris-flow event with no meteorological analogue. Both
simulations peak on 18 Dec 1979 at 38,587 cfs from an ordinary winter storm.
This year fails the timing screen by 152 days and is never adjusted. It is
flagged explicitly in the output so the reason is recorded rather than inferred.

Whether WY1980 belongs in a flood frequency analysis at all is a separate
question this script does not decide.

OUTPUTS
-------
  adjusted_peaks.csv        every shared year, all three peaks, the adjustment,
                            screening results and the reason for each decision
  adjusted_peaks_ssp.csv    WY and adjusted peak only, for HEC-SSP import
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
# Years failing a screen: False = carry through unadjusted, True = omit entirely
DROP_FAILED_SCREEN = False
# Years whose USGS peak is below this are reported separately -- low-flow years
# are where spurious timing mismatches cluster, because the annual maximum is
# not a distinct storm.
LOW_PEAK_CFS = 20000.0

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
    colors = ["#16a085" if p else "#c0392b" for p in table["screen_passed"]]
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
                 "describe the same storm?")
    handles = [Line2D([], [], color="#16a085", lw=6, label="Same event, eligible"),
               Line2D([], [], color="#c0392b", lw=6, label="Different event, not adjusted")]
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
        if PEAK_METHOD == "event":
            p_wcm, t_wcm = event_peak(wcm, t_usgs, EVENT_WINDOW_DAYS)
            p_obs, t_obs = event_peak(obs, t_usgs, EVENT_WINDOW_DAYS)
            if not (np.isfinite(p_wcm) and np.isfinite(p_obs)):
                rows.append({
                    "WY": wy, "usgs": p_usgs, "wcm": np.nan, "obs": np.nan,
                    "delta_wcm_minus_obs": np.nan, "adjusted_peak": p_usgs,
                    "adjusted": False, "screen_passed": False,
                    "spread_days": np.nan, "t_usgs": t_usgs.date(),
                    "t_wcm": pd.NaT, "t_obs": pd.NaT,
                    "low_peak_year": p_usgs < LOW_PEAK_CFS,
                    "decision": ("no adjustment: no simulated data within %d days "
                                 "of the observed peak %s -- the Obs_RC window "
                                 "does not cover this storm"
                                 % (EVENT_WINDOW_DAYS, t_usgs.date()))})
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
        adjust = bool(passed and delta > 0)
        if adjust:
            decision = "adjusted +%.0f cfs" % delta
        elif passed:
            decision = "no adjustment: Obs_RC peak >= WCM_RC peak"
        else:
            decision = "no adjustment: %s" % reason
        rows.append({
            "WY": wy,
            "usgs": p_usgs, "wcm": p_wcm, "obs": p_obs,
            "delta_wcm_minus_obs": delta,
            "adjusted_peak": p_usgs + delta if adjust else p_usgs,
            "adjusted": adjust,
            "screen_passed": passed,
            "spread_days": spread,
            "t_usgs": t_usgs.date(), "t_wcm": t_wcm.date(), "t_obs": t_obs.date(),
            "low_peak_year": p_usgs < LOW_PEAK_CFS,
            "decision": decision,
        })

    table = pd.DataFrame(rows)
    if DROP_FAILED_SCREEN:
        table = table[table["screen_passed"]].reset_index(drop=True)

    table.to_csv(OUT_CSV, index=False, float_format="%.1f")
    table[["WY", "adjusted_peak"]].to_csv(OUT_SSP_CSV, index=False,
                                          float_format="%.0f")

    if os.path.exists(OUT_DSS):
        os.remove(OUT_DSS)
    dss = HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION)
    try:
        span = pd.date_range(pd.Timestamp(int(table["WY"].min()) - 1, 10, 1),
                             pd.Timestamp(int(table["WY"].max()), 9, 30), freq="D")
        values = pd.Series(SENTINEL, index=span)
        for _, row in table.iterrows():
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
    print("NOT ADJUSTED, Obs_RC peak >= WCM_RC peak : %d years" % len(no_gain))
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
    print("-" * 78)
    print("Peaks CSV   : %s" % OUT_CSV)
    print("SSP CSV     : %s" % OUT_SSP_CSV)
    print("DSS         : %s  %s" % (OUT_DSS, OUT_DSS_PATH))
    print("Plots       : %s.png and %s_screening.png" % (PLOT_STEM, PLOT_STEM))


main()
