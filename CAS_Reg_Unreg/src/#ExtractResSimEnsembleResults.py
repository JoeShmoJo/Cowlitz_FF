"""
extract_results.py
"""

import csv
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ---------------------------------------------------------------------------
# USER CONFIG
# ---------------------------------------------------------------------------
# EXTERNAL: requires the ResSim simulation output (not in this repository)
SIM_DSS     = r"C:/Projects/Cowlitz_Flow_Frequency/ResSim/NWP_CowlitzLewis/watershed/NWP_CowlitzLewis/rss/Unreg_2009_2025/simulation.dss"
import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

OBS_DSS     = r"../../CAS_Unreg_FF/data/obsData.dss"
OUT_DSS     = r"../data/results.dss"
MAPPING_CSV = r"../output/ensemble_mapping.csv"

OBS_START = "01Oct2009"
OBS_END   = "30sep2025"

SENTINEL     = -901.0
SENTINEL_TOL = 0.5
# First day of Water year 2010

# Provide one example member pathname per output you want.
# The suffix after | (ENSEMBLE--0 or ENSEMBLE--1) is parsed automatically
# and used to filter so the correct simulation is used.
SIM_RESULTS_PATH_DICT = [
    ("//CASTLEROCK_NWS/FLOW/01JAN2000/1HOUR/C:000001|ENSEMBLE--1/",        "ROUTE_TEST_RESULTS"),
    ("//CASTLEROCK_NWS/FLOW-UNREG/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",  "UNREG_RESULTS"),
    ("//CASTLEROCK_NWS/FLOW/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",        "REG_RESULTS"),
    ("//MAYFIELD_OUT/FLOW-UNREG/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",    "UNREG_RESULTS"),
    ("//MOSSYROCK_OUT/FLOW-UNREG/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",   "UNREG_RESULTS"),
    ("//MAYFIELD-POOL/ELEV/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",         "REG_RESULTS"),
    ("//MOSSYROCK-POOL/ELEV/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",        "REG_RESULTS"),
    ("//MAYFIELD-POOL/FLOW-OUT/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",     "REG_RESULTS"),
    ("//MOSSYROCK-POOL/FLOW-OUT/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",    "REG_RESULTS"),
    ("//MAYFIELD-POOL/FLOW-IN/01JAN2000/1HOUR/C:000001|ENSEMBLE--0/",      "REG_RESULTS"),
]

OBS_PATH_DICT = [
    ("/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW/*/1HOUR/USGS/", "OBS_CASTLEROCK_FLOW"),
    ("//MOS/FLOW-OUT/*/1HOUR/CWMS/",                                   "OBS_MOS_FLOW"),
    ("//MOS/ELEV/*/1HOUR/CWMS/",                                       "OBS_MOS_ELEV"),
    ("//MAY/FLOW-OUT_PEAKCLEAN_2009_2026/*/1HOUR/CWMS/",               "OBS_MAY_FLOW"),
    ("//MAY/ELEV/*/1HOUR/CWMS/",                                       "OBS_MAY_ELEV"),
    ("//MAY/FLOW-IN-CALC-RAW-PEAKS/*/1HOUR/CWMS/",                     "OBS_MAY_FLOW"),
    ("//MOS/ELEV-RULECURVE/*/1HOUR/CENWP-CALC/",                       "OBS_MOS_RULECURVE"),
]

UNITS_BY_PART_C = {
    "ELEV":       "FEET",
    "FLOW-OUT":   "CFS",
    "FLOW-UNREG": "CFS",
    "FLOW":       "CFS",
}

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def parse_pathname(pathname):
    parts = pathname.split("/")
    if len(parts) >= 7:
        return parts[1:7]  # [A, B, C, D, E, F]
    return None

def extract_member_number(part_f):
    """Extract integer member number from e.g. 'C:000001|ENSEMBLE--0' -> 1"""
    if part_f.startswith("C:") and "|" in part_f:
        try:
            return int(part_f[2:part_f.index("|")])
        except ValueError:
            pass
    return None

def extract_ens_suffix(part_f):
    """Extract suffix after | e.g. 'C:000001|ENSEMBLE--0' -> 'ENSEMBLE--0'"""
    if "|" in part_f:
        return part_f[part_f.index("|") + 1:]
    return None

def fmt_dt(dt):
    try:
        return dt.strftime("%-d%b%Y %H%M").upper()
    except ValueError:
        return dt.strftime("%#d%b%Y %H%M").upper()

def fmt_irr_dt(dt):
    """Format datetime as pydsstools irregular time string: '02JUL2010 1200'"""
    return dt.strftime("%d%b%Y %H%M").upper()

def make_part_d(start, end):
    return f"{fmt_dt(start)} - {fmt_dt(end)}"

def load_mapping(csv_path):
    mapping = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            num   = int(row["member_number"])
            start = datetime.strptime(row["actual_start"], "%Y-%m-%d %H:%M")
            end   = datetime.strptime(row["actual_end"],   "%Y-%m-%d %H:%M")
            mapping[num] = (start, end)
    return mapping

def safe_read_ts(src_file, pathname):
    tsc  = src_file.read_ts(pathname.strip())
    vals = np.atleast_1d(np.array(tsc.values, dtype=float))
    if vals.size == 0:
        raise ValueError("Empty array returned")
    return vals, tsc

def water_year(dt):
    return dt.year + 1 if dt.month >= 10 else dt.year

# ---------------------------------------------------------------------------
# WATER YEAR PEAKS
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# WATER YEAR PEAKS
# ---------------------------------------------------------------------------

def compute_wy_peaks(src_part_b, src_part_c, part_f_src):
    """
    Read the full regular 1HOUR record using a wildcard D-part,
    compute WY max 1/24/72/120-hour averages, and write IR-Year records.
    """

    print("\nComputing water-year peaks...")
    print("-" * 60)

    dss = HecDss.Open(OUT_DSS, version=6)

    try:
        src_path = f"//{src_part_b}/{src_part_c}/*/1HOUR/{part_f_src}/"
        print(f"  Reading: {src_path}")

        src_tsc = dss.read_ts(src_path)

        vals = np.atleast_1d(np.array(src_tsc.values, dtype=float))
        units = src_tsc.units or "CFS"
        part_a = src_tsc.pathname.split("/")[1] if src_tsc.pathname else ""

        interval_raw = int(src_tsc.interval)

        if interval_raw == 60:
            interval_hours = 1.0        # minutes
        elif interval_raw == 3600:
            interval_hours = 1.0        # seconds
        else:
            raise ValueError(f"Expected 1HOUR interval, got src_tsc.interval={interval_raw}")
        
        start_dt = None
        for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M"):
            try:
                start_dt = datetime.strptime(src_tsc.startDateTime.title(), fmt)
                break
            except ValueError:
                continue

        if start_dt is None:
            raise ValueError(f"Cannot parse startDateTime: {src_tsc.startDateTime!r}")

        bad = (
            np.isclose(vals, SENTINEL, atol=SENTINEL_TOL)
            | (vals < 0.0)
            | (~np.isfinite(vals))
        )

        vals_fixed = np.where(bad, -10000.0, vals.astype(float))

        idx = pd.date_range(
            start=start_dt,
            periods=len(vals_fixed),
            freq="1h"
        )

        series = pd.Series(vals_fixed, index=idx)

        duration_config = [
            ("1HR", 1),
            ("24HR", 24),
            ("72HR", 72),
            ("120HR", 120),
        ]

        for label, window in duration_config:

            if window == 1:
                rolled = series.copy()
            else:
                rolled = series.rolling(
                    window=window,
                    min_periods=window
                ).mean()

            rolled = rolled.dropna()

            wy_peaks = {}

            for ts, val in rolled.items():

                if window == 1:
                    peak_time = ts.to_pydatetime()
                else:
                    window_start = ts.to_pydatetime() - timedelta(hours=window - 1)
                    window_end = ts.to_pydatetime()
                    peak_time = window_start + (window_end - window_start) / 2

                wy = water_year(peak_time)

                if wy not in wy_peaks or float(val) > wy_peaks[wy][1]:
                    wy_peaks[wy] = (peak_time, float(val))

            wy_peaks = {
                wy: tv
                for wy, tv in wy_peaks.items()
                if tv[1] > -9999.0
            }

            if not wy_peaks:
                print(f"  {label}: no valid peaks found")
                continue

            sorted_wys = sorted(wy_peaks)
            peak_times = [wy_peaks[w][0] for w in sorted_wys]
            peak_vals = [wy_peaks[w][1] for w in sorted_wys]

            out_part_d = f"{fmt_dt(peak_times[0])} - {fmt_dt(peak_times[-1])}"

            out_path = (
                f"/{part_a}/{src_part_b}/{src_part_c}"
                f"/{out_part_d}/IR-Year/UNREG_RESULTS - WY {label}/"
            )

            out_tsc = TimeSeriesContainer()
            out_tsc.pathname = out_path
            out_tsc.interval = -1
            out_tsc.times = peak_times
            out_tsc.values = peak_vals
            out_tsc.numberValues = len(peak_vals)
            out_tsc.units = units
            out_tsc.type = "INST-VAL"

            dss.deletePathname(out_path)
            dss.put_ts(out_tsc)

            print(f"  {label}: {len(peak_vals)} WY peaks  ->  {out_path}")
            for w, (t, v) in sorted(wy_peaks.items()):
                print(f"    WY{w}: {fmt_irr_dt(t)}  {v:>10.1f} {units}")

    finally:
        dss.close()

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print(f"Simulation DSS   : {SIM_DSS}")
    print(f"Observations DSS : {OBS_DSS}")
    print(f"Output DSS       : {OUT_DSS}")
    print(f"Mapping CSV      : {MAPPING_CSV}")
    print("=" * 70 + "\n")

    mapping = load_mapping(MAPPING_CSV)
    print(f"Loaded {len(mapping)} member mappings.\n")

    # ── 1. Copy observation paths ────────────────────────────────────────────
    print("Copying observation paths...")
    print("-" * 60)
    with HecDss.Open(OBS_DSS) as obs, HecDss.Open(OUT_DSS) as dst:
        for pattern, out_label in OBS_PATH_DICT:
            print(f"  {pattern}  ->  F={out_label}")
            try:
                matches = obs.getPathnameList(pattern, sort=1)
                if not matches:
                    print(f"    WARNING: no paths found\n")
                    continue

                actual_path = matches[0]
                print(f"    Found {len(matches)} blocks, sample: {actual_path}")

                src_parts = actual_path.split("/")
                # src_parts = ['', A, B, C, D, E, F, '']
                read_path = (f"/{src_parts[1]}/{src_parts[2]}/{src_parts[3]}"
                             f"/{OBS_START} - {OBS_END}"
                             f"/{src_parts[5]}/{src_parts[6]}/")

                tsc  = obs.read_ts(read_path.strip())
                vals = np.atleast_1d(np.array(tsc.values, dtype=float))
                n_valid = int(np.sum(~np.isclose(vals, SENTINEL, atol=SENTINEL_TOL)))

                out_path = (f"/{src_parts[1]}/{src_parts[2]}/{src_parts[3]}"
                            f"/{OBS_START} - {OBS_END}"
                            f"/{src_parts[5]}/{out_label}/")

                out_tsc = TimeSeriesContainer()
                out_tsc.pathname      = out_path
                out_tsc.startDateTime = tsc.startDateTime
                out_tsc.numberValues  = len(vals)
                out_tsc.interval      = tsc.interval // 60
                out_tsc.values        = vals.tolist()
                out_tsc.units         = tsc.units
                out_tsc.type          = tsc.type

                dst.put_ts(out_tsc)
                print(f"    n={len(vals)}  valid={n_valid}  -> {out_path}")
                print(f"    Written OK.\n")
            except Exception as e:
                import traceback
                print(f"    ERROR: {e}")
                traceback.print_exc()
                print()

    # ── 2. Extract and combine simulation results ────────────────────────────
    print("\nExtracting simulation results...")
    print("-" * 60)
    with HecDss.Open(SIM_DSS) as sim, HecDss.Open(OUT_DSS) as dst:
        for example_path, part_f_out in SIM_RESULTS_PATH_DICT:
            parts = parse_pathname(example_path)
            if parts is None:
                print(f"  ERROR: could not parse: {example_path}\n")
                continue
            part_a, part_b, part_c, _, part_e, example_f = parts

            # Extract ensemble suffix from the example Part F
            ens_suffix = extract_ens_suffix(example_f)
            if ens_suffix is None:
                print(f"  ERROR: could not extract suffix from F={example_f}\n")
                continue

            print(f"\nProcessing: /{part_a}/{part_b}/{part_c}/  "
                  f"suffix={ens_suffix}  ->  F={part_f_out}")
            print("-" * 60)

            pat = f"/{part_a}/{part_b}/{part_c}/*/{part_e}/*/"
            try:
                all_paths = sim.getPathnameList(pat, sort=1)
            except Exception as e:
                print(f"  ERROR getting path list: {e}")
                continue

            if not all_paths:
                print(f"  WARNING: no paths matched {pat}")
                continue

            # Build member number -> pathname lookup filtered by ens_suffix
            path_by_member = {}
            for p in all_paths:
                p_parts = parse_pathname(p)
                if p_parts is None:
                    continue
                part_f = p_parts[5]
                # Only include paths matching the correct ensemble suffix
                if extract_ens_suffix(part_f) != ens_suffix:
                    continue
                num = extract_member_number(part_f)
                if num is not None:
                    path_by_member[num] = p

            print(f"  Members found ({ens_suffix}): {sorted(path_by_member.keys())}")

            combined = {}
            units    = UNITS_BY_PART_C.get(part_c, "CFS")

            for member_num in sorted(mapping.keys()):
                actual_start, actual_end = mapping[member_num]

                if member_num not in path_by_member:
                    print(f"  Member {member_num:>3}: not in sim output (skipping)")
                    continue

                try:
                    vals, tsc = safe_read_ts(sim, path_by_member[member_num])
                    if tsc.units:
                        units = tsc.units
                except Exception as e:
                    print(f"  Member {member_num:>3}: READ ERROR -> {e}")
                    continue

                n_valid  = int(np.sum(~np.isclose(vals, SENTINEL, atol=SENTINEL_TOL)))
                inserted = 0
                for i, v in enumerate(vals):
                    if not np.isclose(v, SENTINEL, atol=SENTINEL_TOL):
                        ts = actual_start + timedelta(hours=i)
                        combined[ts] = float(v)
                        inserted += 1

                print(f"  Member {member_num:>3}: {actual_start} -> {actual_end}  "
                      f"valid={n_valid}/{len(vals)}  inserted={inserted}")

            if not combined:
                print(f"  ERROR: no valid data collected.")
                continue

            all_ts    = sorted(combined.keys())
            out_start = all_ts[0]
            out_end   = all_ts[-1]
            n_out     = int((out_end - out_start).total_seconds() / 3600) + 1

            out_vals = np.full(n_out, SENTINEL)
            for ts, v in combined.items():
                idx = int((ts - out_start).total_seconds() / 3600)
                out_vals[idx] = v

            n_filled   = int(np.sum(~np.isclose(out_vals, SENTINEL, atol=SENTINEL_TOL)))
            out_part_d = make_part_d(out_start, out_end)
            out_path   = f"/{part_a}/{part_b}/{part_c}/{out_part_d}/{part_e}/{part_f_out}/"

            print(f"\n  Combined record:")
            print(f"    {out_start} -> {out_end}  "
                  f"({n_out} slots, {n_filled} filled)")
            print(f"    units={units}")
            print(f"    -> {out_path}")

            out_tsc = TimeSeriesContainer()
            out_tsc.pathname      = out_path
            out_tsc.startDateTime = out_start.strftime("%d%b%Y %H:%M:%S").upper()
            out_tsc.numberValues  = n_out
            out_tsc.interval      = 60
            out_tsc.values        = out_vals.tolist()
            out_tsc.units         = units
            out_tsc.type          = "INST-VAL"

            dst.put_ts(out_tsc)
            print(f"    Written OK.")

    # ── 3. Water-year peaks from CASTLEROCK_NWS FLOW-UNREG ──────────────────
    compute_wy_peaks(
        src_part_b = "CASTLEROCK_NWS",
        src_part_c = "FLOW-UNREG",
        part_f_src = "UNREG_RESULTS",
    )

    print("\n" + "=" * 70)
    print(f"Done. Results written to {OUT_DSS}")


if __name__ == "__main__":
    main()