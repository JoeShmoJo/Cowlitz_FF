#Extract_Ensemble_To_Timeseries.py
# -*- coding: utf-8 -*-
"""
Reassemble ResSim ensemble output into real-dated continuous time series.

ResSim writes one record per ensemble member, all stamped on the same synthetic
calendar window with an F-part like C:000007|ENSEMBLE--0. This script reads the
members listed in the mapping CSV, moves each one back to the real dates its
member represents, and writes one continuous 1HOUR record per location.

It is deliberately generic: point SIM_DSS and MAPPING_CSV at whichever run you
want, list the locations in RECORDS, and set ENS_SUFFIX to match the alternative
(ENSEMBLE--0, ENSEMBLE--1, ...). Nothing else is ensemble-specific.

Members must not overlap in real time. The 01 Oct -> 01 May windows do not, but
the script checks and reports collisions rather than silently overwriting.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
# EXTERNAL: the ResSim simulation output (too large for the repository)
SIM_DSS = r"C:/Projects/2026_Cowlitz_Flow_Frequency/ResSim/NWP_CowlitzLewis/watershed/NWP_CowlitzLewis/rss/1999.10.02-1200/simulation.dss"
SIM_DSS_VERSION = 6

MAPPING_CSV = r"../output/ensemble_unreg_mapping.csv"
OUT_DSS = r"../output/ResSim_Unreg_Results.dss"
OUT_DSS_VERSION = 6
SUMMARY_CSV = r"../output/diagnostics/ensemble_results_summary.csv"

# F-part suffix after the pipe. ENSEMBLE--0 is the first alternative in the run.
ENS_SUFFIX = "ENSEMBLE--0"

# Locations to pull back. (part_a, part_b, part_c, units, out_f_part)
RECORDS = [
    ("", "CASTLEROCK_NWS",  "FLOW",     "CFS",  "UNREG_RESSIM"),
    ("", "MOSSYROCK-POOL",  "FLOW-OUT", "CFS",  "UNREG_RESSIM"),
]

SENTINEL = -901.0
SENTINEL_TOL = 0.5
DROP_NEGATIVE = False    # True = treat negative values as missing

# ----------------------------------------------------------------------------


def first_stamp(ts):
    """First timestamp of a DSS series, across pydsstools versions."""
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M", "%d%b%Y %H%M%S", "%d%b%Y %H%M",
                "%d %B %Y %H:%M:%S", "%d %B %Y %H:%M"):
        try:
            return pd.Timestamp(datetime.strptime(text, fmt))
        except ValueError:
            continue
    return pd.Timestamp(text)


def series_step(ts, pathname):
    """Time step of a DSS regular series. ts.interval is in seconds."""
    seconds = int(getattr(ts, "interval", 0) or 0)
    if seconds > 0:
        return pd.Timedelta(seconds=seconds)
    e_part = pathname.split("/")[5].upper()
    lookup = {"1MIN": "1min", "15MIN": "15min", "30MIN": "30min",
              "1HOUR": "1h", "6HOUR": "6h", "12HOUR": "12h", "1DAY": "1D"}
    return pd.Timedelta(lookup.get(e_part, "1h"))


def read_member(dss, part_a, part_b, part_c, member):
    """Read one ensemble member as a Series on the SYNTHETIC calendar."""
    f_part = "C:%06d|%s" % (member, ENS_SUFFIX)
    pathname = "/%s/%s/%s/*/1HOUR/%s/" % (part_a, part_b, part_c, f_part)
    ts = dss.read_ts(pathname)
    values = np.atleast_1d(np.array(ts.values, dtype=float))
    nodata = np.atleast_1d(np.array(ts.nodata, dtype=bool))
    values[nodata] = np.nan
    values[np.isclose(values, SENTINEL, atol=SENTINEL_TOL)] = np.nan
    values[values <= -900.0] = np.nan
    if DROP_NEGATIVE:
        values[values < 0.0] = np.nan
    step = series_step(ts, pathname)
    index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    return pd.Series(values, index=index)


def load_mapping(csv_path):
    """Mapping rows written by #Create_Unreg_Ensembles.py."""
    table = pd.read_csv(csv_path, parse_dates=["real_start", "real_end",
                                               "ensemble_start"])
    return table.sort_values("member").reset_index(drop=True)


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


def reassemble(dss, mapping, part_a, part_b, part_c):
    """Move every member back to its real dates; returns a real-dated Series."""
    pieces, misses, collisions = [], [], 0
    for _, row in mapping.iterrows():
        member = int(row["member"])
        try:
            synth = read_member(dss, part_a, part_b, part_c, member)
        except Exception as exc:
            misses.append((member, str(exc).strip()[:60]))
            continue
        ens_start = pd.Timestamp(row["ensemble_start"])
        n_hours = int(row["hours"])
        window = synth.loc[ens_start:ens_start + pd.Timedelta(hours=n_hours - 1)]
        if len(window) < n_hours:
            misses.append((member, "short window: %d of %d" % (len(window), n_hours)))
            continue
        real_index = pd.date_range(pd.Timestamp(row["real_start"]),
                                   periods=n_hours, freq="h")
        pieces.append(pd.Series(window.values[:n_hours], index=real_index))

    if not pieces:
        return None, misses, 0
    combined = pd.concat(pieces)
    collisions = int(combined.index.duplicated().sum())
    if collisions:
        combined = combined[~combined.index.duplicated(keep="first")]
    return combined.sort_index(), misses, collisions


def to_continuous(series):
    """Fill the gaps between members so DSS gets one regular record."""
    full = pd.date_range(series.index[0], series.index[-1], freq="h")
    return series.reindex(full)


def main():
    for path in (os.path.dirname(OUT_DSS), os.path.dirname(SUMMARY_CSV)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    mapping = load_mapping(MAPPING_CSV)
    print("=" * 78)
    print("Simulation DSS : %s" % SIM_DSS)
    print("Mapping        : %s  (%d members)" % (MAPPING_CSV, len(mapping)))
    print("F-part suffix  : %s" % ENS_SUFFIX)
    print("Output DSS     : %s" % OUT_DSS)
    print("=" * 78)

    summary = []
    src = HecDss.Open(SIM_DSS, version=SIM_DSS_VERSION)
    dst = HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION)
    try:
        for part_a, part_b, part_c, units, out_f in RECORDS:
            print("\n%s/%s" % (part_b, part_c))
            series, misses, collisions = reassemble(src, mapping, part_a,
                                                    part_b, part_c)
            if series is None:
                print("   no members read -- check the pathname parts and ENS_SUFFIX")
                for member, why in misses[:5]:
                    print("     member %d: %s" % (member, why))
                continue

            n_members = len(mapping) - len(misses)
            full = to_continuous(series)
            out_values = np.where(np.isfinite(full.values), full.values, SENTINEL)
            # DSS stamps are end-of-period; first value's stamp is one step later
            start_time = (full.index[0] + pd.Timedelta(hours=1)).strftime(
                "%d%b%Y %H:%M:%S").upper()
            pathname = "/%s/%s/%s//1HOUR/%s/" % (part_a, part_b, part_c, out_f)
            dst.put_ts(build_container(pathname, out_values, start_time,
                                       units, "INST-VAL", 60))

            print("   members read : %d of %d" % (n_members, len(mapping)))
            print("   real span    : %s -> %s" % (full.index[0].date(),
                                                  full.index[-1].date()))
            print("   values       : %d slots, %d valid, %d gap"
                  % (len(full), int(np.isfinite(full.values).sum()),
                     int((~np.isfinite(full.values)).sum())))
            print("   range        : %.1f to %.1f %s"
                  % (np.nanmin(full.values), np.nanmax(full.values), units))
            print("   written      : %s" % pathname)
            if collisions:
                print("   WARNING: %d overlapping timestamps, first value kept"
                      % collisions)
            for member, why in misses[:5]:
                print("   MISSING member %d: %s" % (member, why))
            if len(misses) > 5:
                print("   ... and %d more missing members" % (len(misses) - 5))

            summary.append({"part_b": part_b, "part_c": part_c,
                            "members_read": n_members,
                            "members_missing": len(misses),
                            "start": full.index[0], "end": full.index[-1],
                            "valid_hours": int(np.isfinite(full.values).sum()),
                            "gap_hours": int((~np.isfinite(full.values)).sum()),
                            "min": float(np.nanmin(full.values)),
                            "max": float(np.nanmax(full.values)),
                            "overlaps": collisions,
                            "out_pathname": pathname})
    finally:
        src.close()
        dst.close()

    pd.DataFrame(summary).to_csv(SUMMARY_CSV, index=False)
    print("\n" + "-" * 78)
    print("Summary CSV: %s" % SUMMARY_CSV)


main()
