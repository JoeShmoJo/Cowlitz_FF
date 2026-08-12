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
# --- which ensemble set is being reassembled ---------------------------------
# ResSim_WCM_RC  : WCM rule curve run, 01 Oct -> 01 May windows
# ResSim_Obs_RC  : observed rule curve run, 31-day windows on the rising limb
SET_NAME = "ResSim_WCM_RC"  # one of the keys in CONFIG_BY_SET below

# Everything that differs between the runs lives in ONE table. SIM_DSS belongs
# here too: pairing one run's simulation with the other run's mapping produces a
# plausible-looking but completely wrong series, because it stamps run A's
# members onto run B's dates. That is not hypothetical -- it happened.
# EXTERNAL: the ResSim simulation output (too large for the repository)
RSS_ROOT = (r"C:\Projects\2026_Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis"
            r"\watershed\NWP_CowlitzLewis_ResSim4\rss")

CONFIG_BY_SET = {
    "ResSim_WCM_RC": {
        "mapping": r"../output/ensemble_unreg_mapping.csv",
        "sim_dss": RSS_ROOT + r"\WCM_RC\simulation.dss",
    },
    "ResSim_Obs_RC": {
        "mapping": r"../output/ensemble_obs_rc_mapping.csv",
        "sim_dss": RSS_ROOT + r"\OBS_RC\simulation.dss",
    },
    # Synthetics: each member has its own synthetic water year (1801+), so the
    # reassembled series is a run of non-overlapping blocks that the mapping CSV
    # navigates by event / magnitude / pool basis.
    "ResSim_Synth": {
        "mapping": r"../output/ensemble_synthetic_mapping.csv",
        "sim_dss": r"../output/simulation.dss",
    },
}
MAPPING_CSV = CONFIG_BY_SET[SET_NAME]["mapping"]
SIM_DSS = CONFIG_BY_SET[SET_NAME]["sim_dss"]
SIM_DSS_VERSION = None   # None = detect from the file header
OUT_DSS = r"../output/%s.dss" % SET_NAME
OUT_DSS_VERSION = 6
SUMMARY_CSV = r"../output/diagnostics/%s_summary.csv" % SET_NAME

# F-part suffix after the pipe, as ResSim writes it in the simulation output.
# NOTE: the INPUT ensembles are now tagged C:00000N|OBS_RC and C:00000N|WCM_RC,
# but ResSim replaces whatever follows the pipe on OUTPUT with its own
# alternative tag. Check the actual F-part in simulation.dss before running:
# it may be "ENSEMBLE--0", or it may now carry the input tag through.
# ResSim echoes the INPUT records back under the input's own suffix
# (C:000001|SYNTH-Ensemble--0) but writes its COMPUTED results under a plain
# C:000001|Ensemble--0. The computed ones are what get reassembled.
ENS_SUFFIX = "Ensemble--0"

# Locations to pull back. (part_a, part_b, part_c, units, out_f_part)
# B and C parts are matched CASE-INSENSITIVELY against the file's catalog, so
# these do not have to match ResSim's capitalisation.
RECORDS = [
    ("", "CastleRock_NWS", "Flow",       "CFS",  SET_NAME),
    ("", "Mossyrock-Pool", "Flow-OUT",   "CFS",  SET_NAME),
    ("", "Mossyrock-Pool", "Elev",       "FEET", SET_NAME),
    # passed straight through ResSim -- reassemble them to check the mapping
    ("", "Mossyrock-Pool", "Flow-IN",    "CFS",  SET_NAME),
    ("", "CastleRock_NWS", "Flow-Local", "CFS",  SET_NAME),
    ("", "CastleRock_NWS", "Flow-UNREG", "CFS",  SET_NAME),
]

# Round-trip check: reassembled record vs. the record it was built from.
# (out_part_b, out_part_c, source_dss, source_pathname)
CHECK_AGAINST = [] if SET_NAME == "ResSim_Synth" else [
    ("Mossyrock-Pool", "Flow-IN", r"../output/ResSimInflows.dss",
     "//MOSSYROCK/FLOW-IN/*/1HOUR/FOR_RESSIM/"),
    ("CastleRock_NWS", "Flow-Local", r"../output/ResSimInflows.dss",
     "//CASTLE ROCK/FLOW-LOCAL/*/1HOUR/FOR_RESSIM/"),
]
SOURCE_DSS_VERSION = 6
CHECK_TOLERANCE_ABS = 0.5      # cfs; below this is DSS single-precision noise
CHECK_TOLERANCE_REL = 0.0001   # or this fraction of the record peak, whichever is larger

# A member returning more than this multiple of its mapped length means the
# simulation and the ensemble do not belong together.
WINDOW_TOLERANCE = 1.5

SENTINEL = -901.0
SENTINEL_TOL = 0.5
DROP_NEGATIVE = False    # True = treat negative values as missing

# ----------------------------------------------------------------------------


def dss_version(path):
    """DSS file version from the header: byte 12 is 6 for v6, 0 for v7.

    The ResSim runs are a mix -- simulation.dss is v7 while the older ensembles
    are v6 -- and pydsstools needs the version passed explicitly on Linux.
    """
    with open(path, "rb") as handle:
        head = handle.read(16)
    if len(head) < 16 or head[:4] != b"ZDSS":
        return None
    return 6 if head[12] == 6 else 7


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


def catalog_paths(dss):
    """Every pathname in the file, across pydsstools versions.

    The catalog API is not stable between builds: some expose search_path,
    some path_dict, older ones getPathnameList. Try each rather than assume --
    this sandbox has search_path while the user's build does not.
    """
    getter = getattr(dss, "search_path", None)
    if getter is not None:
        try:
            return list(getter("/*/*/*/*/*/*/"))
        except Exception:
            pass
    getter = getattr(dss, "path_dict", None)
    if getter is not None:
        try:
            grouped = getter("/*/*/*/*/*/*/")
            paths = []
            for value in grouped.values():
                paths.extend(list(value))
            if paths:
                return paths
        except Exception:
            pass
    for name in ("getPathnameList", "getCatalogedPathnames", "get_pathnames"):
        getter = getattr(dss, name, None)
        if getter is None:
            continue
        for args in (("/*/*/*/*/*/*/",), ()):
            try:
                return list(getter(*args))
            except Exception:
                continue
    raise RuntimeError(
        "Could not list pathnames: this pydsstools build exposes none of "
        "search_path / path_dict / getPathnameList. Available: %s"
        % ", ".join(m for m in dir(dss) if "path" in m.lower()))


def build_catalog(dss):
    """Index the file's real pathnames by (B, C, member) folded to lower case.

    DSS lookups here are CASE SENSITIVE and ResSim does not preserve the case
    used on input -- a run may hold /CastleRock_NWS/Flow/ while an earlier one
    held /CASTLEROCK_NWS/FLOW/. Constructing pathnames therefore fails silently
    on a case change. Resolving them from the catalog does not.
    """
    index = {}
    suffixes = {}
    for path in catalog_paths(dss):
        parts = path.split("/")
        if len(parts) < 8:
            continue
        f_part = parts[6]
        if "|" not in f_part or not f_part.upper().startswith("C:"):
            continue
        member_text, suffix = f_part.split("|", 1)
        try:
            member = int(member_text[2:])
        except ValueError:
            continue
        key = (parts[2].strip().lower(), parts[3].strip().lower(), member, suffix.lower())
        index.setdefault(key, path)
        suffixes[suffix] = suffixes.get(suffix, 0) + 1
    return index, suffixes


def resolve_suffix(suffixes, wanted):
    """Pick the F-part suffix to use, preferring an exact match on ENS_SUFFIX.

    ResSim carries the INPUT tag through onto the records it echoes back (here
    C:000001|SYNTH-Ensemble--0) but writes its own COMPUTED results under a
    different suffix (C:000001|Ensemble--0). Taking the most common suffix would
    pick the wrong family, so prefer the configured one, then the one attached
    to the most records that is not an echo of the input tag.
    """
    if wanted and wanted.lower() in {k.lower() for k in suffixes}:
        for key in suffixes:
            if key.lower() == wanted.lower():
                return key
    return max(suffixes, key=lambda k: suffixes[k]) if suffixes else wanted


def read_member(dss, catalog, suffix, part_b, part_c, member):
    """Read one ensemble member as a Series on the SYNTHETIC calendar."""
    pathname = catalog.get((part_b.strip().lower(), part_c.strip().lower(),
                            member, suffix.lower()))
    if pathname is None:
        raise KeyError("no /%s/%s/ for member %d with suffix %s"
                       % (part_b, part_c, member, suffix))
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


def reassemble(dss, catalog, suffix, mapping, part_a, part_b, part_c):
    """Move every member back to its real dates; returns a real-dated Series.

    ResSim decides its own simulation window, which is usually NOT the window the
    ensemble was written on -- a run folder named 1999.10.02-1200 produces output
    stamped 02Oct1999 - 29Apr2000 even though the input covered 01Oct - 01May. So
    do not slice to the input window. Instead preserve each value's OFFSET from
    ensemble_start and add it to real_start, which works for any run window:

        real_time = real_start + (output_time - ensemble_start)
    """
    pieces, misses, offsets = [], [], []
    for _, row in mapping.iterrows():
        member = int(row["member"])
        try:
            synth = read_member(dss, catalog, suffix, part_b, part_c, member)
        except Exception as exc:
            misses.append((member, str(exc).strip()[:60]))
            continue
        synth = synth.dropna()          # drops DSS block padding as well as gaps
        if synth.empty:
            misses.append((member, "no valid values"))
            continue
        ens_start = pd.Timestamp(row["ensemble_start"])
        real_start = pd.Timestamp(row["real_start"])
        lead = synth.index[0] - ens_start
        offsets.append((member, lead, synth.index[-1] - ens_start, len(synth)))
        pieces.append(pd.Series(synth.values, index=real_start + (synth.index - ens_start)))

    if not pieces:
        return None, misses, 0, offsets
    combined = pd.concat(pieces)
    collisions = int(combined.index.duplicated().sum())
    if collisions:
        combined = combined[~combined.index.duplicated(keep="first")]
    return combined.sort_index(), misses, collisions, offsets


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
    version = SIM_DSS_VERSION or dss_version(SIM_DSS)
    src = HecDss.Open(SIM_DSS, version=version)
    catalog, suffixes = build_catalog(src)
    suffix = resolve_suffix(suffixes, ENS_SUFFIX)
    print("F-part suffixes in the file: %s"
          % ", ".join("%s (%d)" % (k, v) for k, v in
                      sorted(suffixes.items(), key=lambda kv: -kv[1])[:6]))
    print("Using suffix   : %s%s"
          % (suffix, "" if suffix.lower() == (ENS_SUFFIX or "").lower()
             else "   (ENS_SUFFIX was %r -- not present, fell back)" % ENS_SUFFIX))
    dst = HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION)
    try:
        for part_a, part_b, part_c, units, out_f in RECORDS:
            print("\n%s/%s" % (part_b, part_c))
            series, misses, collisions, offsets = reassemble(
                src, catalog, suffix, mapping, part_a, part_b, part_c)
            if series is None:
                print("   no members read -- check the pathname parts and ENS_SUFFIX")
                for member, why in misses[:5]:
                    print("     member %d: %s" % (member, why))
                continue

            leads = set(o[1] for o in offsets)
            tails = set(o[2] for o in offsets)
            member_hours = int(mapping["hours"].iloc[0])
            got_hours = max(o[3] for o in offsets)
            if got_hours > member_hours * WINDOW_TOLERANCE:
                print("   *** STOP: members returned up to %d hours but the mapping"
                      " says %d ***" % (got_hours, member_hours))
                print("   The simulation window does not match this ensemble.")
                print("   Usual cause: SIM_DSS is a different run than SET_NAME=%s"
                      % SET_NAME)
                print("   SIM_DSS = %s" % SIM_DSS)
                print("   Skipping this record -- the result would be wrong, not"
                      " merely short.")
                continue
            print("   run window   : starts %s after ensemble_start, ends %s after"
                  % (sorted(leads)[0], sorted(tails)[-1]))
            if len(leads) > 1 or len(tails) > 1:
                print("   NOTE: members do not share one run window "
                      "(%d distinct starts, %d distinct ends)" % (len(leads), len(tails)))

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

    if CHECK_AGAINST:
        run_checks()

    print("\n" + "-" * 78)
    print("Summary CSV: %s" % SUMMARY_CSV)


def read_whole(dss_file, pathname, version):
    """Read a full record as a real-dated Series (used for the round-trip check)."""
    dss = HecDss.Open(dss_file, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.atleast_1d(np.array(ts.values, dtype=float))
        nodata = np.atleast_1d(np.array(ts.nodata, dtype=bool))
        values[nodata] = np.nan
        values[np.isclose(values, SENTINEL, atol=SENTINEL_TOL)] = np.nan
        values[values <= -900.0] = np.nan
        step = series_step(ts, pathname)
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index()


def run_checks():
    """Compare pass-through records against the source they were built from."""
    print("\n" + "=" * 78)
    print("ROUND-TRIP CHECK -- reassembled vs. original")
    print("=" * 78)
    for out_b, out_c, src_dss, src_path in CHECK_AGAINST:
        out_path = None
        for part_a, part_b, part_c, units, out_f in RECORDS:
            if part_b == out_b and part_c == out_c:
                out_path = "/%s/%s/%s/*/1HOUR/%s/" % (part_a, part_b, part_c, out_f)
        if out_path is None:
            print("\n%s/%s: not in RECORDS, skipping" % (out_b, out_c))
            continue
        try:
            got = read_whole(OUT_DSS, out_path, OUT_DSS_VERSION).dropna()
            want = read_whole(src_dss, src_path, SOURCE_DSS_VERSION).dropna()
        except Exception as exc:
            print("\n%s/%s: could not read -- %s" % (out_b, out_c, exc))
            continue
        shared = got.index.intersection(want.index)
        print("\n%s/%s" % (out_b, out_c))
        print("   reassembled : %d values  %s -> %s"
              % (len(got), got.index[0].date(), got.index[-1].date()))
        print("   original    : %d values  %s -> %s"
              % (len(want), want.index[0].date(), want.index[-1].date()))
        print("   overlapping : %d hours" % len(shared))
        if len(shared) == 0:
            print("   *** NO OVERLAP -- the mapping is wrong ***")
            continue
        diff = (got.loc[shared] - want.loc[shared]).abs()
        scale = max(float(want.loc[shared].abs().max()), 1.0)
        tol = max(CHECK_TOLERANCE_ABS, CHECK_TOLERANCE_REL * scale)
        n_over = int((diff > tol).sum())
        print("   max abs diff: %.4f   mean abs diff: %.6f   (%.4f%% of peak)"
              % (diff.max(), diff.mean(), 100.0 * diff.max() / scale))
        print("   values off by more than %.3f: %d of %d (%.3f%%)"
              % (tol, n_over, len(shared), 100.0 * n_over / len(shared)))
        if n_over == 0:
            print("   OK -- differences are storage precision, timing is aligned")
        else:
            worst = diff.idxmax()
            print("   CHECK: worst at %s -- reassembled %.2f vs original %.2f"
                  % (worst, got.loc[worst], want.loc[worst]))
            off = diff[diff > tol]
            print("   first 3 offenders: %s"
                  % ", ".join("%s (%.2f)" % (t.strftime("%Y-%m-%d %H:%M"), d)
                              for t, d in off.head(3).items()))
            print("   A handful of large diffs usually means a value edit "
                  "(clipping, cleaning) between the source and the ensemble;")
            print("   a systematic offset across many hours means the timing "
                  "mapping is wrong.")


main()
