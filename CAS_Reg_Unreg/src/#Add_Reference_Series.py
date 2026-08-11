#Add_Reference_Series.py
# -*- coding: utf-8 -*-
"""
Add the WCM rule curve and the observed pool to the reassembled result files.

Both are INPUTS, not ResSim output, so there is no reason to pull them back out
of simulation.dss -- where the ensemble F-part is written inconsistently
(C:000001| in some records, C:000053|-ENSEMBLE--0 in others). Instead:

    rule curve       generated from RULE_CURVE_ANCHORS, which is a formula
    observed pool    read from //MOS/ELEV//1DAY/USGS/ and resampled to hourly

Both are then masked to the member windows in the mapping CSV, so they only
exist where simulated results exist, and written into the result DSS with the
same F-part as everything else in that file.

Run this AFTER #Extract_Ensemble_To_Timeseries.py for the set in question.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from datetime import datetime
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
OBS_DSS = r"../../CAS_Unreg_FF/data/obsData.dss"
OBS_DSS_VERSION = 6
PATH_MOS_ELEV_DAILY = "//MOS/ELEV/*/1DAY/USGS/"

# Which result files get what. rule_curve is written to both sets; observed pool
# only to the set that actually starts from it.
TARGETS = [
    {"set_name": "ResSim_WCM_RC",
     "result_dss": r"../output/ResSim_WCM_RC.dss",
     "mapping": r"../output/ensemble_unreg_mapping.csv",
     "rule_curve": True, "observed_pool": False},
    {"set_name": "ResSim_Obs_RC",
     "result_dss": r"../output/ResSim_Obs_RC.dss",
     "mapping": r"../output/ensemble_obs_rc_mapping.csv",
     "rule_curve": True, "observed_pool": True},
]
RESULT_DSS_VERSION = 6

# Same anchors as both ensemble builders: (month, day, elevation)
RULE_CURVE_ANCHORS = [
    (1, 1, 745.5),
    (1, 31, 745.5),
    (6, 1, 778.5),
    (9, 30, 778.5),
    (12, 1, 745.5),
]
RULE_CURVE_PART_B = "MOS"
RULE_CURVE_PART_C = "ELEV-RULECURVE"
OBS_POOL_PART_B = "MOS"
OBS_POOL_PART_C = "ELEV-OBS"

ELEV_DAILY_TO_HOURLY = "interpolate"   # "interpolate" or "step"
# True  = only where simulated results exist (member windows)
# False = continuous over the whole span, gaps and all
MASK_TO_MEMBER_WINDOWS = True
PAD_HOURS = 24            # extend each member window this much on each side

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
    return pd.Series(values, index=index).sort_index()


def daily_to_hourly(daily, how):
    """Daily elevation to hourly. Daily values are labelled on their own day."""
    hourly_index = pd.date_range(daily.index[0],
                                 daily.index[-1] + pd.Timedelta(hours=23), freq="h")
    if how == "step":
        return daily.reindex(hourly_index, method="ffill")
    anchored = daily.copy()
    anchored.index = anchored.index + pd.Timedelta(hours=12)
    return (anchored.reindex(anchored.index.union(hourly_index))
            .interpolate(method="time").reindex(hourly_index))


def rule_curve_on_index(index):
    """Seasonal WCM rule curve for a DatetimeIndex, linear between anchors.

    Interpolate on float HOURS from a fixed epoch. Do NOT mix Index.asi8 with
    Timestamp.value -- pandas 3 returns microseconds from asi8 and nanoseconds
    from .value, and the 1000x mismatch silently flattens the curve.
    """
    epoch = pd.Timestamp(1900, 1, 1)
    knots = []
    for year in range(int(index.year.min()) - 1, int(index.year.max()) + 2):
        for month, day, elev in RULE_CURVE_ANCHORS:
            hours = (pd.Timestamp(year, month, day) - epoch) / pd.Timedelta(hours=1)
            knots.append((hours, float(elev)))
    knots.sort()
    kt = np.array([k[0] for k in knots], dtype=float)
    kv = np.array([k[1] for k in knots], dtype=float)
    x = (index - epoch) / pd.Timedelta(hours=1)
    return np.interp(np.asarray(x, dtype=float), kt, kv)


def member_mask(index, mapping):
    """True where the hour falls inside a member window (padded)."""
    mask = np.zeros(len(index), dtype=bool)
    pad = pd.Timedelta(hours=PAD_HOURS)
    for _, row in mapping.iterrows():
        a = pd.Timestamp(row["real_start"]) - pad
        b = pd.Timestamp(row["real_end"]) + pad
        mask |= (index >= a) & (index <= b)
    return mask


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


def write_series(dss, part_b, part_c, series, units, f_part):
    """Write one continuous 1HOUR record, sentinel in the gaps."""
    full = series.reindex(pd.date_range(series.index[0], series.index[-1], freq="h"))
    values = np.where(np.isfinite(full.values), full.values, SENTINEL)
    start_time = (full.index[0] + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()
    pathname = "//%s/%s//1HOUR/%s/" % (part_b, part_c, f_part)
    dss.put_ts(build_container(pathname, values, start_time, units, "INST-VAL", 60))
    return pathname, full, int(np.isfinite(full.values).sum())


def main():
    obs_daily = read_dss_series(OBS_DSS, PATH_MOS_ELEV_DAILY, OBS_DSS_VERSION).dropna()
    obs_hourly = daily_to_hourly(obs_daily, ELEV_DAILY_TO_HOURLY)
    print("Observed pool record: %s .. %s" % (obs_daily.index[0].date(),
                                              obs_daily.index[-1].date()))

    for target in TARGETS:
        result_dss = target["result_dss"]
        if not os.path.isfile(result_dss):
            print("\n%s: %s not found -- run the extractor for this set first"
                  % (target["set_name"], result_dss))
            continue
        mapping = pd.read_csv(target["mapping"],
                              parse_dates=["real_start", "real_end"])
        span_a = pd.Timestamp(mapping["real_start"].min()) - pd.Timedelta(hours=PAD_HOURS)
        span_b = pd.Timestamp(mapping["real_end"].max()) + pd.Timedelta(hours=PAD_HOURS)
        index = pd.date_range(span_a, span_b, freq="h")
        mask = member_mask(index, mapping) if MASK_TO_MEMBER_WINDOWS else np.ones(len(index), bool)

        print("\n== %s ==" % target["set_name"])
        print("   %d members, span %s .. %s, %d hours inside member windows"
              % (len(mapping), span_a.date(), span_b.date(), int(mask.sum())))

        dss = HecDss.Open(result_dss, version=RESULT_DSS_VERSION)
        try:
            if target["rule_curve"]:
                rc = pd.Series(rule_curve_on_index(index), index=index)
                rc[~mask] = np.nan
                path, full, n = write_series(dss, RULE_CURVE_PART_B, RULE_CURVE_PART_C,
                                             rc.dropna(), "FEET", target["set_name"])
                print("   rule curve   %-46s %d values, %.2f to %.2f ft"
                      % (path, n, np.nanmin(full.values[np.isfinite(full.values)]),
                         np.nanmax(full.values[np.isfinite(full.values)])))
            if target["observed_pool"]:
                pool = obs_hourly.reindex(index)
                pool[~mask] = np.nan
                if pool.notna().sum() == 0:
                    print("   observed pool: no overlap with the member windows")
                else:
                    path, full, n = write_series(dss, OBS_POOL_PART_B, OBS_POOL_PART_C,
                                                 pool.dropna(), "FEET", target["set_name"])
                    print("   observed pool %-46s %d values, %.2f to %.2f ft"
                          % (path, n, np.nanmin(full.values[np.isfinite(full.values)]),
                             np.nanmax(full.values[np.isfinite(full.values)])))
        finally:
            dss.close()

    print("\nDone. Both series are generated from inputs, so this can be re-run")
    print("any time without touching simulation.dss.")


main()
