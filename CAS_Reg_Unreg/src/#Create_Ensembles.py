"""
create_ensemble_dss.py
"""

import sys
import csv
import numpy as np
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

OBS_DSS     = r"../../CAS_Unreg_FF/data/obsData.dss"
# EXTERNAL: requires the ResSim watershed (not in this repository)
OUT_DSS     = r"C:\Projects\Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis\watershed\NWP_CowlitzLewis\shared\ensemble.dss"
MAPPING_CSV = r"../output/ensemble_mapping.csv"

ENS_START    = datetime(2000, 1, 1, 12, 0)   # hour-beginning convention
ENS_END      = datetime(2000, 1, 31, 12, 0)
WINDOW_HOURS = int((ENS_END - ENS_START).total_seconds() / 3600) + 1  # 721

SENTINEL     = -901.0
SENTINEL_TOL = 0.5

VOLCOR_PATH = "//MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/"

UNITS_BY_PART_C = {
    "ELEV":                            "FEET",
    "FLOW-LOCAL-SHAPED":               "CFS",
    "FLOW-IN-CALC-CLEANED-VOLCOR":     "CFS",
    "FLOW-LOCAL":                      "CFS",
    "FLOW":                            "CFS",
    "FLOW-OUT_PEAKCLEAN_2009_2026":    "CFS",
}

SOURCE_PATHS = [
    ("",           "MAY",       "ELEV",                          "1HOUR", "//MAY/ELEV//1HOUR/CWMS/"),
    ("",           "MOS",       "ELEV",                          "1HOUR", "//MOS/ELEV//1HOUR/CWMS/"),
    ("",           "MAY",       "FLOW-LOCAL-SHAPED",             "1HOUR", "//MAY/FLOW-LOCAL-SHAPED//1HOUR/CWMS/"),
    ("",           "MOS",       "FLOW-IN-CALC-CLEANED-VOLCOR",   "1HOUR", "//MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/"),
    ("",           "MAY",       "FLOW-OUT_PEAKCLEAN_2009_2026",  "1HOUR", "//MAY/FLOW-OUT_PEAKCLEAN_2009_2026//1HOUR/CWMS/"),
    ("",           "MOS",       "FLOW-OUT",  "1HOUR", "//MOS/FLOW-OUT//1HOUR/CWMS/"),
    ("COWLITZ RIVER AT CASTLE ROCK, WA", "14243000", "FLOW-LOCAL","1HOUR", "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW-LOCAL//1HOUR/USGS/"),
    ("ZERO_DUMMY", "ZERO",       "FLOW",                          "1HOUR", "/ZERO_DUMMY/ZERO/FLOW/01SEP2008 - 01SEP2026/1HOUR/DUMMY/"),
    ("TOUTLE RIVER AT TOWER ROAD NEAR SILVER LAKE, WA", "14242580", "FLOW", "1HOUR", "/TOUTLE RIVER AT TOWER ROAD NEAR SILVER LAKE, WA/14242580/FLOW//1HOUR/USGS/"),
    ("COWLITZ RIVER AT CASTLE ROCK, WA", "14243000", "FLOW",      "1HOUR", "/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW//1HOUR/USGS/"),
]

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def parse_dss_datetime(s):
    """
    pydsstools returns startDateTime as hour-ENDING of the first value.
    Subtract 1 hour to get the hour-BEGINNING label we use internally.
    e.g. '01Sep2008 01:00:00' means value covers 00:00-01:00 -> label as 00:00
    """
    s = s.strip()
    for fmt in ["%d%b%Y %H:%M:%S", "%d%b%Y %H:%M", "%d%b%Y %H%M%S", "%d%b%Y %H%M"]:
        try:
            return datetime.strptime(s, fmt) - timedelta(hours=1)
        except ValueError:
            continue
    raise ValueError(f"Cannot parse startDateTime: {s!r}")

def fmt_dt(dt):
    try:
        return dt.strftime("%-d%b%Y %H%M").upper()
    except ValueError:
        return dt.strftime("%#d%b%Y %H%M").upper()

def ensemble_part_f(n):
    return f"C:{n:06d}|"

def ensemble_part_d(start, end):
    return f"{fmt_dt(start)} - {fmt_dt(end)}"

def is_sentinel(v):
    return abs(float(v) - SENTINEL) <= SENTINEL_TOL

def read_raw(src_file, path):
    tsc = src_file.read_ts(path.strip())
    vals = np.array(tsc.values, dtype=float)
    start_dt = parse_dss_datetime(tsc.startDateTime)   # now hour-beginning
    interval_h = tsc.interval / 3600.0
    print(f"    raw={tsc.startDateTime!r}  parsed(hr-beg)={start_dt}  "
          f"interval_h={interval_h}")
    return vals, start_dt, interval_h

def make_datetime_index(start_dt, n, interval_h):
    step = timedelta(hours=interval_h)
    return [start_dt + i * step for i in range(n)]

def find_period_starts(vals, times):
    starts = []
    in_period = False
    for i, v in enumerate(vals):
        if not is_sentinel(v) and not in_period:
            starts.append(times[i])
            in_period = True
        elif is_sentinel(v) and in_period:
            in_period = False
    return starts

def extract_window(vals, times, win_start):
    ts_map = {t: v for t, v in zip(times, vals)}
    step   = timedelta(hours=1)
    out    = np.full(WINDOW_HOURS, np.nan)
    for i in range(WINDOW_HOURS):
        ts = win_start + i * step
        v  = ts_map.get(ts, np.nan)
        if not np.isnan(v) and not is_sentinel(v):
            out[i] = v
    last_valid = np.nan
    for i in range(WINDOW_HOURS):
        if np.isnan(out[i]):
            out[i] = last_valid
        else:
            last_valid = out[i]
    out = np.where(np.isnan(out), SENTINEL, out)
    return out

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------


ens_d = ensemble_part_d(ENS_START, ENS_END)

print("=" * 70)
print(f"Ensemble window : {ENS_START}  ->  {ENS_END}  ({WINDOW_HOURS} slots)")
print(f"Source DSS      : {OBS_DSS}")
print(f"Output DSS      : {OUT_DSS}")
print(f"Mapping CSV     : {MAPPING_CSV}")
print("=" * 70 + "\n")

# ── 1. Read VOLCOR ───────────────────────────────────────────────────────
print(f"Reading VOLCOR:\n  {VOLCOR_PATH}")
with HecDss.Open(OBS_DSS) as src:
    volcor_vals, volcor_start, volcor_ih = read_raw(src, VOLCOR_PATH)

volcor_times  = make_datetime_index(volcor_start, len(volcor_vals), volcor_ih)
period_starts = find_period_starts(volcor_vals, volcor_times)

print(f"\nFound {len(period_starts)} continuous periods:")
for i, ps in enumerate(period_starts, 1):
    ev_end = ps + timedelta(hours=WINDOW_HOURS - 1)
    print(f"  Period {i:>3}: {ps}  ->  {ev_end}")
print()

if not period_starts:
    sys.exit("ERROR: No valid periods found.")

# ── 2. Write mapping CSV ─────────────────────────────────────────────────
with open(MAPPING_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["member_number", "actual_start", "actual_end",
                        "ens_start", "ens_end"])
    for i, ps in enumerate(period_starts, 1):
        ev_end = ps + timedelta(hours=WINDOW_HOURS - 1)
        writer.writerow([
            i,
            ps.strftime("%Y-%m-%d %H:%M"),
            ev_end.strftime("%Y-%m-%d %H:%M"),
            ENS_START.strftime("%Y-%m-%d %H:%M"),
            ENS_END.strftime("%Y-%m-%d %H:%M"),
        ])
print(f"Mapping CSV written: {MAPPING_CSV}\n")

# ── 3. Read all source paths ─────────────────────────────────────────────
print("Reading all source paths into memory...")
all_data = []
with HecDss.Open(OBS_DSS) as src:
    for part_a, part_b, part_c, part_e, full_path in SOURCE_PATHS:
        print(f"  {full_path}")
        try:
            vals, start_dt, ih = read_raw(src, full_path)
            times   = make_datetime_index(start_dt, len(vals), ih)
            n_valid = sum(1 for v in vals if not is_sentinel(v))
            first_valid_idx = next(
                (i for i, v in enumerate(vals) if not is_sentinel(v)), None)
            if first_valid_idx is not None:
                print(f"    First valid: idx={first_valid_idx}  "
                        f"time={times[first_valid_idx]}  "
                        f"val={vals[first_valid_idx]:.3f}")
            print(f"    OK  n={len(vals)}  valid={n_valid}  "
                    f"start={times[0]}  end={times[-1]}")
            all_data.append((part_a, part_b, part_c, part_e, vals, times))
        except Exception as e:
            print(f"    WARNING: could not read -> {e}")
            all_data.append((part_a, part_b, part_c, part_e, None, None))

# Cross-check at Period 1 start
print(f"\nCross-check at Period 1 start ({period_starts[0]}):")
for part_a, part_b, part_c, part_e, vals, times in all_data:
    if vals is None:
        continue
    ts_map = {t: v for t, v in zip(times, vals)}
    v = ts_map.get(period_starts[0], "NOT IN INDEX")
    print(f"  {part_b}/{part_c}: {v}")
print()

# ── 4. Write ensemble members ────────────────────────────────────────────
print(f"Writing ensemble members to {OUT_DSS}\n")

# ENS_START is hour-beginning; DSS write expects hour-ending -> add 1 hour
dss_write_start = (ENS_START + timedelta(hours=1)).strftime(
                        "%d%b%Y %H:%M:%S").upper()
print(f"  DSS startDateTime for write: {dss_write_start}\n")

with HecDss.Open(OUT_DSS, version = 6) as dst:
    for member_num, ps in enumerate(period_starts, start=1):
        ev_end = ps + timedelta(hours=WINDOW_HOURS - 1)
        pf     = ensemble_part_f(member_num)
        print(f"  Member {member_num:>3}: {ps} -> {ev_end}")

        for part_a, part_b, part_c, part_e, vals, times in all_data:
            if vals is None:
                print(f"    SKIP: {part_b}/{part_c}")
                continue

            arr     = extract_window(vals, times, ps)
            n_valid = int(np.sum(arr != SENTINEL))

            if member_num == 1:
                print(f"    {part_b}/{part_c}: "
                        f"[0]={arr[0]:.2f}  [1]={arr[1]:.2f}  [2]={arr[2]:.2f}")

            pathname = f"/{part_a}/{part_b}/{part_c}/{ens_d}/{part_e}/{pf}/"
            units    = UNITS_BY_PART_C.get(part_c, "CFS")

            tsc = TimeSeriesContainer()
            tsc.pathname      = pathname
            tsc.startDateTime = dss_write_start   # hour-ending for DSS
            tsc.numberValues  = WINDOW_HOURS
            tsc.interval      = 60
            tsc.values        = arr.tolist()
            tsc.units         = units
            tsc.type          = "INST-VAL"

            dst.put_ts(tsc)
        print()

n_good = sum(1 for d in all_data if d[4] is not None)
print("=" * 70)
print(f"Done. {len(period_starts)} members x {n_good} paths = "
        f"{len(period_starts) * n_good} records written.")
print(f"Mapping saved to: {MAPPING_CSV}")




# ---------------------------------------------------------------------------
# DIAGNOSTIC CSVs — timestamp alignment inspection
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# DIAGNOSTIC CSVs — timestamp alignment inspection
# ---------------------------------------------------------------------------
import csv as csv_module

print("\nWriting diagnostic CSVs...")

# Which all_data indices to inspect
INSPECT = [
    ("CASTLEROCK_FLOW",  8),   # /COWLITZ.../14243000/FLOW/
    ("MAY_PEAKCLEAN",    4),   # //MAY/FLOW-OUT_PEAKCLEAN/
    ("VOLCOR",          -1),   # special: volcor
]

# Write one CSV per inspected path
for label, idx in INSPECT:
    if idx == -1:
        vals_i  = volcor_vals
        times_i = volcor_times
        pb, pc  = "MOS", "FLOW-IN-CALC-CLEANED-VOLCOR"
    else:
        part_a, part_b, part_c, part_e, vals_i, times_i = all_data[idx]
        pb, pc = part_b, part_c

    if vals_i is None:
        print(f"  SKIP {label} — no data")
        continue

    out_csv = rf"../output/diagnostics/inspect_{label}.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv_module.writer(f)
        writer.writerow(["slot_index", "datetime", "value", "note"])

        # First 30 slots
        for i in range(min(30, len(times_i))):
            writer.writerow([i, times_i[i], vals_i[i], "START"])

        # Last 30 slots
        for i in range(max(0, len(times_i) - 30), len(times_i)):
            writer.writerow([i, times_i[i], vals_i[i], "END"])

        # For each period: 10 slots before and 20 slots after period start
        for member_num, ps in enumerate(period_starts, start=1):
            # Find exact index of period start
            exact_idx = None
            for i, t in enumerate(times_i):
                if t == ps:
                    exact_idx = i
                    break

            if exact_idx is None:
                # Find nearest instead
                diffs = [abs((t - ps).total_seconds()) for t in times_i]
                exact_idx = int(np.argmin(diffs))
                nearest   = times_i[exact_idx]
                offset_h  = (nearest - ps).total_seconds() / 3600
                note_pfx  = f"PERIOD_{member_num}_NEAREST(off={offset_h:+.1f}h)"
            else:
                note_pfx = f"PERIOD_{member_num}_EXACT"

            for i in range(max(0, exact_idx - 10), min(len(times_i), exact_idx + 20)):
                marker = " <-- PERIOD START" if i == exact_idx else ""
                writer.writerow([i, times_i[i], vals_i[i],
                                 f"{note_pfx}{marker}"])

    print(f"  Written: {out_csv}")

# Master alignment CSV — for every period, show what value each path
# has at the period start timestamp, and what slot that maps to in each path
align_csv = rf"../output/diagnostics/alignment_check.csv"
with open(align_csv, "w", newline="") as f:
    writer = csv_module.writer(f)

    # Header
    header = ["member", "period_start"]
    for part_a, part_b, part_c, part_e, vals_i, times_i in all_data:
        header += [f"{part_b}/{part_c}_idx",
                   f"{part_b}/{part_c}_time",
                   f"{part_b}/{part_c}_val"]
    writer.writerow(header)

    for member_num, ps in enumerate(period_starts, start=1):
        row = [member_num, ps]
        for part_a, part_b, part_c, part_e, vals_i, times_i in all_data:
            if vals_i is None:
                row += ["NO DATA", "NO DATA", "NO DATA"]
                continue
            # Find index of ps in this path's time index
            exact_idx = None
            for i, t in enumerate(times_i):
                if t == ps:
                    exact_idx = i
                    break
            if exact_idx is None:
                diffs     = [abs((t - ps).total_seconds()) for t in times_i]
                exact_idx = int(np.argmin(diffs))
                nearest   = times_i[exact_idx]
                offset_h  = (nearest - ps).total_seconds() / 3600
                row += [f"NEAREST(idx={exact_idx},off={offset_h:+.1f}h)",
                        nearest, vals_i[exact_idx]]
            else:
                row += [exact_idx, times_i[exact_idx], vals_i[exact_idx]]
        writer.writerow(row)

print(f"  Written: {align_csv}")

# Ensemble slot 0 verification CSV
# For member 1: show what value is in slot 0 of the extracted array
# vs what's in the actual DSS file at ENS_START
ens_verify_csv = rf"../output/diagnostics/ensemble_slot0_check.csv"
with open(ens_verify_csv, "w", newline="") as f:
    writer = csv_module.writer(f)
    writer.writerow(["path", "period_start", "slot0_value",
                     "ens_start_written", "note"])

    ps = period_starts[0]
    for part_a, part_b, part_c, part_e, vals_i, times_i in all_data:
        if vals_i is None:
            continue
        arr = extract_window(vals_i, times_i, ps)
        writer.writerow([
            f"{part_b}/{part_c}",
            ps,
            arr[0],
            ENS_START,
            f"slot0 should appear at {ENS_START} in DSS viewer"
        ])
        # Also write slots 0-5 so we can check the sequence
        for i in range(6):
            writer.writerow([
                f"{part_b}/{part_c}",
                ps + timedelta(hours=i),
                arr[i],
                ENS_START + timedelta(hours=i),
                f"slot{i}"
            ])

print(f"  Written: {ens_verify_csv}")
print("\nOpen these CSVs and compare against the DSS viewer to identify the shift.")