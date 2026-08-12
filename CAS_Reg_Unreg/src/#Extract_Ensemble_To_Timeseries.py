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
        "mapping": r"../output/ensemble_wcm_rc_mapping.csv",
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

# ...and a member returning less than this fraction of its mapped length means
# the read is truncated. If EVERY member is short by the same factor the cause
# is structural (a block-level pathname, a wrong mapping), not missing data.
WINDOW_SHORT_FRACTION = 0.5

# --- regulated vs unregulated sanity check -----------------------------------
# The reservoir cannot make a flood bigger. At a flood event the regulated flow
# at Castle Rock must sit at or below the unregulated flow for the same hour,
# and the regulated PEAK must sit below the unregulated peak for the year. The
# legitimate exceptions are narrow: evacuating storage on a falling limb, and
# low-flow periods where a minimum release or a refill drawdown puts more water
# in the river than nature would. Anything else -- above all a regulated peak
# above the unregulated peak at a large event -- means the two records do not
# belong together: a mismatched member, a mis-mapped window, or an operation
# set releasing more than it should.
#
# Checked pairs are (part_b, regulated part_c, unregulated part_c). Both parts
# must appear in RECORDS above; a pair naming a record that was not built is
# skipped with a note rather than failing.
REG_UNREG_CHECKS = [
    ("CastleRock_NWS", "Flow", "Flow-UNREG"),
]
# Ignore differences below this -- DSS is single precision and the two records
# take slightly different routes through ResSim.
REG_UNREG_TOL_CFS = 1.0
# Hourly exceedances where the UNREGULATED flow is below this are expected
# behaviour (minimum release, refill drawdown) and are counted separately rather
# than flagged. Annual PEAK exceedances are flagged regardless of magnitude.
REG_UNREG_LOW_FLOW_CFS = 20000.0
# Contiguous exceedance runs shorter than this are timing noise, not a real
# crossing, and are left out of the episode table.
REG_UNREG_MIN_EPISODE_HOURS = 3
REG_UNREG_WY_CSV = r"../output/diagnostics/%s_reg_vs_unreg_wy.csv" % SET_NAME
REG_UNREG_EPISODE_CSV = (r"../output/diagnostics/%s_reg_gt_unreg_episodes.csv"
                         % SET_NAME)

WATER_YEAR_START_MONTH = 10

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
        # The catalog lists ONE PATHNAME PER STORAGE BLOCK -- a 1HOUR record
        # covering 01 Oct -> 01 May appears eight times, once per month, each
        # with its own D part. Reading a block-specific pathname returns only
        # that block, so blank the D part here and let DSS assemble the whole
        # record. Keeping the first D part seen instead silently returns ~1
        # month per member, and because catalog order differs per record the
        # surviving month differs too -- which is how Flow and Flow-UNREG ended
        # up on non-overlapping dates.
        parts[4] = ""
        path = "/".join(parts)
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
    """Mapping rows written by the matching #Create_*_Ensembles.py."""
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


def water_year(stamp):
    return stamp.year + (1 if stamp.month >= WATER_YEAR_START_MONTH else 0)


def exceedance_episodes(diff, reg, unreg, min_hours):
    """Contiguous hourly runs where regulated exceeds unregulated.

    Works on the gap-filled hourly index, so a NaN hour breaks an episode rather
    than being bridged -- two crossings either side of a data gap are two
    episodes, not one long one.
    """
    over = (diff > REG_UNREG_TOL_CFS).values
    stamps = diff.index
    d_values = diff.values
    r_values = reg.values
    u_values = unreg.values
    episodes = []
    i, n = 0, len(over)
    while i < n:
        if not over[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and over[j + 1]:
            j += 1
        hours = j - i + 1
        if hours >= min_hours:
            k = i + int(np.nanargmax(d_values[i:j + 1]))
            episodes.append({
                "WY": water_year(stamps[i]),
                "start": stamps[i],
                "end": stamps[j],
                "hours": hours,
                "max_exceed_cfs": float(d_values[k]),
                "max_exceed_time": stamps[k],
                "reg_at_max": float(r_values[k]),
                "unreg_at_max": float(u_values[k]),
                "unreg_below_low_flow": bool(u_values[k] < REG_UNREG_LOW_FLOW_CFS),
            })
        i = j + 1
    return episodes


def check_reg_vs_unreg(built):
    """Flag every hour and every water year where regulated exceeds unregulated.

    Two separate questions, reported separately because they fail for different
    reasons:

      PEAK   the water-year regulated maximum against the water-year
             unregulated maximum. This is the one that matters -- a regulated
             peak above the unregulated peak at a flood event is not physical.
      HOURLY every hour where the regulated flow exceeds the unregulated flow.
             Some of these are legitimate (storage evacuation, minimum release),
             so they are split by whether the unregulated flow at that hour is
             above REG_UNREG_LOW_FLOW_CFS.
    """
    print("\n" + "=" * 78)
    print("REGULATED vs UNREGULATED CHECK")
    print("=" * 78)

    all_rows, all_episodes = [], []
    for part_b, reg_c, unreg_c in REG_UNREG_CHECKS:
        reg_raw = built.get((part_b, reg_c))
        unreg_raw = built.get((part_b, unreg_c))
        if reg_raw is None or unreg_raw is None:
            missing = reg_c if reg_raw is None else unreg_c
            print("\n%s: /%s/ was not built -- add it to RECORDS to run this "
                  "check" % (part_b, missing))
            continue

        span = pd.date_range(min(reg_raw.index[0], unreg_raw.index[0]),
                             max(reg_raw.index[-1], unreg_raw.index[-1]), freq="h")
        reg = reg_raw.reindex(span)
        unreg = unreg_raw.reindex(span)
        diff = reg - unreg                      # NaN where either is missing
        keys = np.array([water_year(t) for t in span])

        print("\n%s:  %s vs %s" % (part_b, reg_c, unreg_c))
        print("   comparable hours: %d of %d"
              % (int(diff.notna().sum()), len(span)))

        rows = []
        for wy in np.unique(keys):
            block = keys == wy
            r_wy = reg[block].dropna()
            u_wy = unreg[block].dropna()
            d_wy = diff[block].dropna()
            if r_wy.empty or u_wy.empty:
                continue
            reg_peak = float(r_wy.max())
            unreg_peak = float(u_wy.max())
            t_reg = r_wy.idxmax()
            t_unreg = u_wy.idxmax()
            peak_diff = reg_peak - unreg_peak
            over = d_wy[d_wy > REG_UNREG_TOL_CFS]
            over_high = over[unreg.reindex(over.index) >= REG_UNREG_LOW_FLOW_CFS]
            rows.append({
                "part_b": part_b,
                "WY": int(wy),
                "reg_peak": reg_peak,
                "unreg_peak": unreg_peak,
                "reg_peak_time": t_reg,
                "unreg_peak_time": t_unreg,
                "peak_offset_hrs": (t_reg - t_unreg).total_seconds() / 3600.0,
                "reg_at_unreg_peak": float(reg.get(t_unreg, np.nan)),
                "peak_diff_cfs": peak_diff,
                "peak_ratio": reg_peak / unreg_peak if unreg_peak > 0 else np.nan,
                "REG_PEAK_EXCEEDS": bool(peak_diff > REG_UNREG_TOL_CFS),
                "hours_reg_gt_unreg": int(len(over)),
                "hours_reg_gt_unreg_above_lowflow": int(len(over_high)),
                "max_hourly_exceed_cfs": float(over.max()) if len(over) else 0.0,
                "max_hourly_exceed_time": over.idxmax() if len(over) else pd.NaT,
                "valid_hours": int(len(d_wy)),
            })

        table = pd.DataFrame(rows)
        if table.empty:
            print("   no overlapping water years")
            continue
        episodes = exceedance_episodes(diff, reg, unreg,
                                       REG_UNREG_MIN_EPISODE_HOURS)
        for ep in episodes:
            ep["part_b"] = part_b
        all_rows.append(table)
        all_episodes.extend(episodes)

        bad = table[table["REG_PEAK_EXCEEDS"]]
        print("   water years compared        : %d" % len(table))
        print("   REG PEAK ABOVE UNREG PEAK   : %d" % len(bad))
        for _, row in bad.sort_values("peak_diff_cfs", ascending=False).iterrows():
            note = "  [low flow year]" if row["unreg_peak"] < REG_UNREG_LOW_FLOW_CFS else ""
            print("      WY%d  reg %.0f > unreg %.0f  (+%.0f cfs, %.2fx)  "
                  "reg %s / unreg %s%s"
                  % (row["WY"], row["reg_peak"], row["unreg_peak"],
                     row["peak_diff_cfs"], row["peak_ratio"],
                     pd.Timestamp(row["reg_peak_time"]).strftime("%Y-%m-%d %H:%M"),
                     pd.Timestamp(row["unreg_peak_time"]).strftime("%Y-%m-%d %H:%M"),
                     note))
        big = bad[bad["unreg_peak"] >= REG_UNREG_LOW_FLOW_CFS]
        if len(big):
            print("      *** %d of these are at unreg peaks >= %.0f cfs -- not a"
                  " low-flow artefact ***" % (len(big), REG_UNREG_LOW_FLOW_CFS))
            print("      The reservoir cannot raise a flood peak. Check that"
                  " SIM_DSS and MAPPING_CSV are the same run, that the")
            print("      unregulated record was reassembled on the same"
                  " members, and the ResSim operation set for that year.")
        elif len(bad):
            print("      all at unreg peaks below %.0f cfs -- consistent with"
                  " minimum release / refill drawdown"
                  % REG_UNREG_LOW_FLOW_CFS)

        total_hours = int(table["hours_reg_gt_unreg"].sum())
        high_hours = int(table["hours_reg_gt_unreg_above_lowflow"].sum())
        print("   hours reg > unreg           : %d (%d with unreg >= %.0f cfs)"
              % (total_hours, high_hours, REG_UNREG_LOW_FLOW_CFS))
        print("   exceedance episodes >= %dh   : %d"
              % (REG_UNREG_MIN_EPISODE_HOURS, len(episodes)))
        worst = sorted(episodes, key=lambda e: -e["max_exceed_cfs"])[:5]
        for ep in worst:
            print("      WY%d  %s -> %s  %dh  max +%.0f cfs "
                  "(reg %.0f vs unreg %.0f)"
                  % (ep["WY"], ep["start"].strftime("%Y-%m-%d %H:%M"),
                     ep["end"].strftime("%Y-%m-%d %H:%M"), ep["hours"],
                     ep["max_exceed_cfs"], ep["reg_at_max"], ep["unreg_at_max"]))

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(REG_UNREG_WY_CSV, index=False, float_format="%.2f")
        print("\n   water-year table : %s" % REG_UNREG_WY_CSV)
    if all_episodes:
        pd.DataFrame(all_episodes).to_csv(REG_UNREG_EPISODE_CSV, index=False,
                                          float_format="%.2f")
        print("   episode table    : %s" % REG_UNREG_EPISODE_CSV)


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
    built = {}          # (part_b, part_c) -> reassembled hourly series
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
            mapped = dict(zip(mapping["member"].astype(int),
                               mapping["hours"].astype(int)))
            short = [(m, n, mapped.get(m, member_hours)) for m, _, _, n in offsets
                     if n < mapped.get(m, member_hours) * WINDOW_SHORT_FRACTION]
            if short and len(short) == len(offsets):
                print("   *** STOP: every member came back short -- e.g. member"
                      " %d returned %d of %d mapped hours ***"
                      % (short[0][0], short[0][1], short[0][2]))
                print("   All members short by the same factor is structural, not"
                      " missing data: usually the pathname carries a D part, so")
                print("   DSS returns a single storage block instead of the whole"
                      " record. Skipping -- the result would be a wrong slice.")
                continue
            if short:
                print("   NOTE: %d of %d members returned under %.0f%% of their"
                      " mapped hours (first: member %d, %d of %d)"
                      % (len(short), len(offsets), 100 * WINDOW_SHORT_FRACTION,
                         short[0][0], short[0][1], short[0][2]))
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
            built[(part_b, part_c)] = full

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

    if REG_UNREG_CHECKS:
        check_reg_vs_unreg(built)

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
