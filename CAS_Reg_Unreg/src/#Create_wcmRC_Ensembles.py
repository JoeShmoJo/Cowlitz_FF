#Create_Unreg_Ensembles.py
# -*- coding: utf-8 -*-
"""
Build ResSim ensemble members from the unregulated ResSim inflows.

One member per water year, each covering 01 Oct -> 01 May, sliced from the two
records in ResSimInflows.dss:

    //MOSSYROCK/FLOW-IN//1HOUR/FOR_RESSIM/
    //CASTLE ROCK/FLOW-LOCAL//1HOUR/FOR_RESSIM/

A Mossyrock rule curve ensemble is generated alongside them from the seasonal
anchor points in RULE_CURVE_ANCHORS, interpolated linearly between anchors and
evaluated on the real calendar dates of each member (so leap years are handled).

Every member is written on the same synthetic D-part window with an ensemble
F-part (C:000001|, C:000002|, ...), which is what ResSim expects.  Because every
member starts on 01 Oct and runs the same number of hours, mapping ResSim output
back to real dates is a fixed offset per member -- ensemble_unreg_mapping.csv
records the real start and end for each.

Outputs: the ensemble DSS (external, in the ResSim watershed) and the mapping CSV.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import csv
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
IN_DSS = r"../output/ResSimInflows.dss"
IN_DSS_VERSION = 6

# EXTERNAL: requires the ResSim watershed (not in this repository)
OUT_DSS = r"../output/ensemble_wcm_rc.dss"
OUT_DSS_VERSION = 6
MAPPING_CSV = r"../output/ensemble_wcm_rc_mapping.csv"
COVERAGE_CSV = r"../output/diagnostics/ensemble_wcm_rc_coverage.csv"

PATH_MOS_IN = "//MOSSYROCK/FLOW-IN/01JUN1928 - 01MAY2026/1HOUR/FOR_RESSIM/"
PATH_CAS_LOCAL = "//CASTLE ROCK/FLOW-LOCAL/01Jul1928 - 12May2026/1HOUR/FOR_RESSIM/"

# Member window: 01 Oct of the start year through 01 May of the following year.
WINDOW_START_MONTH = 10
WINDOW_START_DAY = 1
WINDOW_END_MONTH = 5
WINDOW_END_DAY = 1

# Synthetic calendar the ensemble members are written on. Any non-leap year
# works; 2000 matches the regulated ensembles already in the watershed.
# Text appended to the F-part after the collection ID, e.g. C:000002|WCM_RC
ENS_SUFFIX = "WCM_RC"

ENS_LABEL_START = datetime(1999, 10, 1, 0, 0)   # hour-beginning

FIRST_YEAR = None        # None = earliest year with full coverage
LAST_YEAR = None         # None = latest year with a complete window
MIN_COVERAGE_PCT = 0.0   # skip members below this % of non-missing hours
CLIP_NEGATIVE_FLOW = False   # MOSSYROCK FLOW-IN contains negative values LEAVE THEM IN - THE CALCULATION IS VOLUME CORRECT

SENTINEL = -901.0
SENTINEL_TOL = 0.5

# Rule curve: (month, day, elevation). Linear between anchors, wrapping Dec->Jan.
RULE_CURVE_ANCHORS = [
    (1, 1, 745.5),
    (1, 31, 745.5),
    (6, 1, 778.5),
    (9, 30, 778.5),
    (12, 1, 745.5),
]
WRITE_RULE_CURVE = True
RULE_CURVE_PART_B = "MOS"
RULE_CURVE_PART_C = "ELEV-RULECURVE"

# (part_a, part_b, part_c, units) for each series written per member
OUT_PARTS_MOS = ("", "MOSSYROCK", "FLOW-IN", "CFS")
OUT_PARTS_CAS = ("", "CASTLE ROCK", "FLOW-LOCAL", "CFS")

# ----------------------------------------------------------------------------


def first_stamp(ts):
    """First timestamp of a DSS series, across pydsstools versions.

    Some builds yield HecTime objects from ts.times, others yield strings, and
    ts.startDateTime is not always populated. Only the first stamp is needed --
    the rest follow from the interval -- so avoid walking the whole generator.
    """
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
    """Read a DSS regular time series into a pandas Series on period-BEGINNING labels."""
    dss = HecDss.Open(dss_file, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        nodata = np.array(ts.nodata, dtype=bool)
        values[nodata] = np.nan
        values[values <= -900.0] = np.nan
        step = series_step(ts, pathname)
        # DSS stamps are END of period; step back one interval to label the start
        start = first_stamp(ts) - step
        index = pd.date_range(start, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index()


def window_hours(start_year):
    """Hour count of the 01 Oct -> 01 May window beginning in start_year, inclusive."""
    a = datetime(start_year, WINDOW_START_MONTH, WINDOW_START_DAY)
    b = datetime(start_year + 1, WINDOW_END_MONTH, WINDOW_END_DAY)
    return int((b - a).total_seconds() // 3600) + 1


def max_window_hours(years):
    """Longest window across the member years, so every member is the same length."""
    return max(window_hours(y) for y in years)


def candidate_years(mos, cas):
    """Years whose full window falls inside both input records."""
    first = max(mos.index[0], cas.index[0])
    last = min(mos.index[-1], cas.index[-1])
    years = []
    for y in range(first.year - 1, last.year + 1):
        a = pd.Timestamp(y, WINDOW_START_MONTH, WINDOW_START_DAY)
        b = a + pd.Timedelta(hours=max_window_hours([y]) - 1)
        if a >= first and b <= last:
            years.append(y)
    if FIRST_YEAR is not None:
        years = [y for y in years if y >= FIRST_YEAR]
    if LAST_YEAR is not None:
        years = [y for y in years if y <= LAST_YEAR]
    return years


def rule_curve_on_index(index):
    """Seasonal rule curve elevation for a DatetimeIndex, linear between anchors.

    Anchors are repeated for every calendar year the index spans (plus one on
    each side) so the Dec -> Jan segment wraps correctly. Interpolating on real
    timestamps means leap years need no special handling.
    """
    # Interpolate on float HOURS from a fixed epoch. Do NOT mix Index.asi8 with
    # Timestamp.value -- pandas 3 returns microseconds from asi8 and nanoseconds
    # from .value, and the 1000x mismatch silently flattens the curve.
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


def slice_member(series, start, n_hours):
    """Values for one member window, carry-forward filled, gaps left as sentinel."""
    idx = pd.date_range(start, periods=n_hours, freq="h")
    vals = series.reindex(idx)
    n_missing = int(vals.isna().sum())
    filled = vals.ffill()
    filled = filled.where(filled.notna(), SENTINEL)
    return filled.values, n_missing, idx


def ensemble_f_part(member):
    return "C:%06d|%s" % (member, ENS_SUFFIX)


def fmt_dss(dt):
    return dt.strftime("%d%b%Y %H%M").upper()


def ensemble_d_part(start, n_hours):
    return "%s - %s" % (fmt_dss(start), fmt_dss(start + timedelta(hours=n_hours - 1)))


def main():
    for path in (os.path.dirname(MAPPING_CSV), os.path.dirname(COVERAGE_CSV)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    mos = read_dss_series(IN_DSS, PATH_MOS_IN, IN_DSS_VERSION)
    cas = read_dss_series(IN_DSS, PATH_CAS_LOCAL, IN_DSS_VERSION)
    if CLIP_NEGATIVE_FLOW:
        n_neg = int((mos < 0).sum() + (cas < 0).sum())
        mos = mos.clip(lower=0.0)
        cas = cas.clip(lower=0.0)
    else:
        n_neg = 0

    years = candidate_years(mos, cas)
    n_hours = max_window_hours(years)
    ens_start = ENS_LABEL_START
    ens_d = ensemble_d_part(ens_start, n_hours)
    # DSS stamps are end-of-period; the first value's stamp is one hour later
    dss_start = (ens_start + timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()

    print("=" * 78)
    print("Window        : %02d/%02d -> %02d/%02d   %d hours (%.1f days) per member"
          % (WINDOW_START_MONTH, WINDOW_START_DAY, WINDOW_END_MONTH, WINDOW_END_DAY,
             n_hours, n_hours / 24.0))
    print("Members       : %d  (%d - %d)" % (len(years), years[0], years[-1]))
    print("F-part        : C:<member>|%s" % ENS_SUFFIX)
    print("Ensemble D-part: %s" % ens_d)
    print("Input  DSS    : %s" % IN_DSS)
    print("Output DSS    : %s" % OUT_DSS)
    if CLIP_NEGATIVE_FLOW and n_neg:
        print("Clipped %d negative flow values to zero" % n_neg)
    print("=" * 78)

    mapping_rows, coverage_rows, written = [], [], 0
    with HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION) as dst:
        member = 0
        for year in years:
            start = pd.Timestamp(year, WINDOW_START_MONTH, WINDOW_START_DAY)
            mos_vals, mos_miss, idx = slice_member(mos, start, n_hours)
            cas_vals, cas_miss, _ = slice_member(cas, start, n_hours)
            cov = 100.0 * (1.0 - max(mos_miss, cas_miss) / float(n_hours))
            coverage_rows.append({"water_year": year + 1, "start": start,
                                  "mos_missing_hours": mos_miss,
                                  "cas_missing_hours": cas_miss,
                                  "coverage_pct": round(cov, 2)})
            if cov < MIN_COVERAGE_PCT:
                print("  skip %d  coverage %.1f%%" % (year, cov))
                continue

            member += 1
            f_part = ensemble_f_part(member)
            real_end = start + pd.Timedelta(hours=n_hours - 1)

            series = [(OUT_PARTS_MOS, mos_vals, "INST-VAL"),
                      (OUT_PARTS_CAS, cas_vals, "INST-VAL")]
            if WRITE_RULE_CURVE:
                rc = rule_curve_on_index(idx)
                series.append((("", RULE_CURVE_PART_B, RULE_CURVE_PART_C, "FEET"),
                               rc, "INST-VAL"))

            for (part_a, part_b, part_c, units), vals, dtype in series:
                pathname = "/%s/%s/%s/%s/1HOUR/%s/" % (part_a, part_b, part_c,
                                                       ens_d, f_part)
                dst.put_ts(build_container(pathname, vals, dss_start, units,
                                           dtype, 60))
                written += 1

            mapping_rows.append({"member": member, "ensemble_f_part": f_part,
                                 "water_year": year + 1,
                                 "real_start": start, "real_end": real_end,
                                 "ensemble_start": ens_start,
                                 "hours": n_hours,
                                 "offset_hours_to_real": 0,
                                 "coverage_pct": round(cov, 2)})
            if member <= 3 or member == len(years):
                print("  member %3d  WY%d  %s -> %s  coverage %.1f%%"
                      % (member, year + 1, start.date(), real_end.date(), cov))

    pd.DataFrame(mapping_rows).to_csv(MAPPING_CSV, index=False)
    pd.DataFrame(coverage_rows).to_csv(COVERAGE_CSV, index=False)

    low = [r for r in coverage_rows if r["coverage_pct"] < 99.9]
    print("-" * 78)
    print("Members written : %d   records written: %d" % (len(mapping_rows), written))
    print("Mapping CSV     : %s" % MAPPING_CSV)
    print("Coverage CSV    : %s" % COVERAGE_CSV)
    print("Members with gaps (carry-forward filled): %d" % len(low))
    for r in low[:30]:
        print("   WY%d  coverage %.1f%%  (MOS %d h, CAS %d h missing)"
              % (r["water_year"], r["coverage_pct"],
                 r["mos_missing_hours"], r["cas_missing_hours"]))
    print("To map ResSim output back: every member starts on its real_start and")
    print("runs %d hours, so slot i of member N is real_start(N) + i hours." % n_hours)


main()
