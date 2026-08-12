#Create_Synthetic_Ensembles.py
# -*- coding: utf-8 -*-
"""
Build scaled synthetic flood ensembles for the unregulated-regulated transform.

WHY
---
Only one water year in 99 exceeds the unregulated 100-year peak, so the top of
the unreg-reg curve is essentially unconstrained by the record. These synthetics
populate that range by scaling observed events up to target magnitudes.

DESIGN -- a full factorial over the three things that could matter:

    SHAPE     4 source events spanning sharp to sustained hydrographs.
              Attenuation among large observed events is NOT predictable from
              shape, magnitude, antecedent flow or timing (all r-squared < 0.10,
              n=18), so shape is varied rather than assumed away.

    MAGNITUDE 4 targets: the 100-year, 250-year and 500-year unregulated peaks,
              plus one BEYOND the 500-year. The last one exists so the drawn
              curve is not extrapolating past its final point exactly where the
              answer is needed.

    POOL      3 starting elevations: the WCM rule curve, the 50% pool duration
              curve for that date, and the observed pool on that date. Written
              as separate members so the results can be compared, or pooled into
              a Monte Carlo over starting conditions.

    -> 4 x 4 x 3 = 48 members.

SCALING
-------
BOTH the Mossyrock inflow and the Castle Rock local are multiplied by the SAME
factor. The local is a median 48% of unregulated flow at large events, and the
project has no control over it, so scaling only the reservoir inflow would
shrink the uncontrolled fraction and change attenuation for reasons unrelated to
magnitude. Uniform scaling preserves the observed coincidence between the two.

The scale factor is set so the scaled unregulated PEAK (Mossyrock inflow plus
local, the same quantity the transform uses on its x-axis) hits the target.

Outputs: the ensemble DSS, a mapping CSV carrying scale factor / target AEP /
source event / pool basis for every member, and a diagnostic plot per event.
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
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
IN_DSS = r"../output/ResSimInflows.dss"
OBS_DSS = r"../../CAS_Unreg_FF/data/obsData.dss"

# EXTERNAL: requires the ResSim watershed (not in this repository)
OUT_DSS = (r"C:\Projects\2026_Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis"
           r"\watershed\NWP_CowlitzLewis\shared\ensemble_synthetic.dss")
OUT_DSS_VERSION = 7

MAPPING_CSV = r"../output/ensemble_synthetic_mapping.csv"
PLOT_STEM = r"../output/diagnostics/ensemble_synthetic"

PATH_MOS_IN = "//MOSSYROCK/FLOW-IN/*/1HOUR/FOR_RESSIM/"
PATH_CAS_LOCAL = "//CASTLE ROCK/FLOW-LOCAL/*/1HOUR/FOR_RESSIM/"
PATH_MOS_ELEV_DAILY = "//MOS/ELEV/*/1DAY/USGS/"

ENS_SUFFIX = "SYNTH"

# --- source events: (label, peak date, shape 5day/peak, note) ---------------
# Window is centred on the peak; WINDOW_BEFORE/AFTER days on each side.
SOURCE_EVENTS = [
    ("Dec1977", "1977-12-02", "sharp     (shape 0.47, rank 2)"),
    ("Feb1996", "1996-02-09", "sharp-mid (shape 0.50, rank 1)"),
    ("Dec1933", "1933-12-22", "sustained (shape 0.70, pass-through outlier)"),
    ("Dec2025", "2025-12-11", "sustained (shape 0.72, recent)"),
]
WINDOW_BEFORE_DAYS = 10       # lead-in, so the pool responds before the peak
WINDOW_AFTER_DAYS = 20        # recession, long enough to pass the flood

# --- target unregulated PEAK magnitudes -------------------------------------
# From CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv, Duration = Peak.
TARGETS = [
    ("100yr", 0.010, 168884.0),
    ("250yr", 0.004, 195000.0),
    ("500yr", 0.002, 211932.0),
    ("beyond", None, 255000.0),   # past the reported curve, to anchor the end slope
]

# --- starting pool bases ----------------------------------------------------
# "rulecurve" : the WCM rule curve on that calendar date
# "duration50": the 50% exceedance pool elevation for that calendar date,
#               computed from the observed daily record
# "observed"  : the observed pool on that calendar date
POOL_BASES = ["rulecurve", "duration50", "observed"]
# The conservative convention from the workflow email is the HIGHER of the
# three. Set this True to add that as a fourth basis.
ADD_HIGHEST_OF_ALL = False

RULE_CURVE_ANCHORS = [
    (1, 1, 745.5), (1, 31, 745.5), (6, 1, 778.5), (9, 30, 778.5), (12, 1, 745.5),
]
LOOKBACK_DAYS = 1
ELEV_DAILY_TO_HOURLY = "interpolate"
DURATION_50_MIN_YEARS = 10     # min years of record behind a duration-curve value

# Synthetic calendar the members are labelled on (any non-leap year)
ENS_LABEL_START = datetime(1999, 10, 1, 0, 0)

# Every member shares its source event's real dates, so the 12 members built
# from one event would collide when results are reassembled. Each member is
# therefore given its own synthetic water year: the month and day of the source
# event are kept, so seasonality stays visible, and only the year is replaced.
# Years start well before the record so they can never be mistaken for real
# ones. The mapping CSV carries the synthetic date AND the true source date.
SYNTH_YEAR_BASE = 1800

CLIP_NEGATIVE_FLOW = True
SENTINEL = -901.0

# ----------------------------------------------------------------------------


def dss_version(path):
    """DSS file version from the header: byte 12 is 6 for v6, 0 for v7."""
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


def read_dss_series(dss_file, pathname):
    """Read a DSS regular time series into a Series on period-BEGINNING labels."""
    version = dss_version(dss_file)
    dss = HecDss.Open(dss_file, version=version) if version else HecDss.Open(dss_file)
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
    """Seasonal WCM rule curve, linear between anchors.

    Interpolate on float HOURS from a fixed epoch. Do NOT mix Index.asi8 with
    Timestamp.value -- pandas 3 returns microseconds from asi8 and nanoseconds
    from .value, and the 1000x mismatch silently flattens the curve.
    """
    epoch = pd.Timestamp(1900, 1, 1)
    knots = []
    for year in range(int(index.year.min()) - 1, int(index.year.max()) + 2):
        for month, day, elev in RULE_CURVE_ANCHORS:
            knots.append(((pd.Timestamp(year, month, day) - epoch)
                          / pd.Timedelta(hours=1), float(elev)))
    knots.sort()
    kt = np.array([k[0] for k in knots], dtype=float)
    kv = np.array([k[1] for k in knots], dtype=float)
    x = (index - epoch) / pd.Timedelta(hours=1)
    return np.interp(np.asarray(x, dtype=float), kt, kv)


def duration_50_by_dayofyear(daily_elev, min_years):
    """50% exceedance pool elevation for each (month, day) from the record."""
    frame = pd.DataFrame({"elev": daily_elev.values,
                          "month": daily_elev.index.month,
                          "day": daily_elev.index.day})
    grouped = frame.groupby(["month", "day"])["elev"]
    median = grouped.median()
    counts = grouped.count()
    return median[counts >= min_years]


def duration_50_on_index(table, index):
    """Look up the 50% duration curve for a DatetimeIndex (Feb 29 -> Feb 28)."""
    out = np.full(len(index), np.nan)
    for i in range(len(index)):
        key = (int(index.month[i]), int(index.day[i]))
        if key in table.index:
            out[i] = table.loc[key]
        elif key == (2, 29) and (2, 28) in table.index:
            out[i] = table.loc[(2, 28)]
    return pd.Series(out, index=index).interpolate(limit_direction="both").values


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


def to_dss_values(values):
    filled = pd.Series(values).ffill().bfill()
    return np.where(np.isfinite(filled.values), filled.values, SENTINEL)


def fmt_dss(stamp):
    return pd.Timestamp(stamp).strftime("%d%b%Y %H%M").upper()


def d_part(start, n_hours):
    return "%s - %s" % (fmt_dss(start),
                        fmt_dss(pd.Timestamp(start) + pd.Timedelta(hours=n_hours - 1)))


def plot_events(events, mapping, stem):
    """Source hydrographs and the scaled family built from each."""
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9), squeeze=False)
    for k, ev in enumerate(events):
        ax = axes[k // 2][k % 2]
        hours = np.arange(len(ev["total"])) / 24.0
        for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
            if row["pool_basis"] != POOL_BASES[0]:
                continue
            ax.plot(hours, ev["total"].values * row["scale_factor"],
                    lw=1.2, label="%s  x%.2f" % (row["target"], row["scale_factor"]))
        ax.plot(hours, ev["total"].values, color="k", lw=1.8, label="observed x1.00")
        ax.plot(hours, ev["cas"].fillna(0.0).values, color="0.6", lw=0.9, ls=":",
                label="local (observed)")
        ax.set_title("%s   %s" % (ev["label"], ev["note"]), fontsize=10)
        ax.set_xlabel("Days into the member window", fontsize=8)
        ax.set_ylabel("Unregulated flow (cfs)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)
    fig.suptitle("Synthetic source events and their scaled family "
                 "(Mossyrock inflow + Castle Rock local)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("%s_events.png" % stem, dpi=150)
    plt.close(fig)


def plot_pools(events, stem):
    """The three starting pool bases for each event date."""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(events))
    width = 0.26
    colors = {"rulecurve": "#8e44ad", "duration50": "#16a085", "observed": "#2c7fb8"}
    for i, basis in enumerate(POOL_BASES):
        vals = [ev["pools"].get(basis, np.nan) for ev in events]
        ax.bar(x + (i - 1) * width, vals, width, color=colors.get(basis, "0.5"),
               label=basis)
    ax.set_xticks(x)
    ax.set_xticklabels([ev["label"] for ev in events])
    ax.set_ylim(700, 790)
    ax.set_ylabel("Starting pool elevation (ft)")
    ax.set_title("Starting pool basis by source event -- each becomes its own "
                 "ensemble member")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig("%s_pools.png" % stem, dpi=150)
    plt.close(fig)


def main():
    for path in (os.path.dirname(MAPPING_CSV), os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    mos = read_dss_series(IN_DSS, PATH_MOS_IN).dropna()
    cas = read_dss_series(IN_DSS, PATH_CAS_LOCAL).dropna()
    if CLIP_NEGATIVE_FLOW:
        mos = mos.clip(lower=0.0)
        cas = cas.clip(lower=0.0)
    elev_daily = read_dss_series(OBS_DSS, PATH_MOS_ELEV_DAILY).dropna()
    elev_hourly = daily_to_hourly(elev_daily, ELEV_DAILY_TO_HOURLY)
    dur50 = duration_50_by_dayofyear(elev_daily, DURATION_50_MIN_YEARS)

    n_hours = (WINDOW_BEFORE_DAYS + WINDOW_AFTER_DAYS) * 24
    elev_hours = n_hours + LOOKBACK_DAYS * 24
    ens_start = pd.Timestamp(ENS_LABEL_START)
    elev_ens_start = ens_start - pd.Timedelta(days=LOOKBACK_DAYS)
    flow_d = d_part(ens_start, n_hours)
    elev_d = d_part(elev_ens_start, elev_hours)
    flow_write = (ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()
    elev_write = (elev_ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()

    bases = list(POOL_BASES) + (["highest"] if ADD_HIGHEST_OF_ALL else [])

    print("=" * 78)
    print("Window     : %d days before / %d days after the peak (%d hours)"
          % (WINDOW_BEFORE_DAYS, WINDOW_AFTER_DAYS, n_hours))
    print("Design     : %d events x %d magnitudes x %d pool bases = %d members"
          % (len(SOURCE_EVENTS), len(TARGETS), len(bases),
             len(SOURCE_EVENTS) * len(TARGETS) * len(bases)))
    print("Scaling    : BOTH Mossyrock inflow and Castle Rock local, same factor")
    print("Output     : %s" % OUT_DSS)
    print("=" * 78)

    events = []
    for label, peak_date, note in SOURCE_EVENTS:
        centre = pd.Timestamp(peak_date)
        start = (centre - pd.Timedelta(days=WINDOW_BEFORE_DAYS)).normalize()
        index = pd.date_range(start, periods=n_hours, freq="h")
        mos_w = mos.reindex(index)
        cas_w = cas.reindex(index)
        total = mos_w.fillna(0.0) + cas_w.fillna(0.0)
        elev_index = pd.date_range(start - pd.Timedelta(days=LOOKBACK_DAYS),
                                   periods=elev_hours, freq="h")
        obs_pool = elev_hourly.reindex(elev_index)
        rc_pool = pd.Series(rule_curve_on_index(elev_index), index=elev_index)
        d50_pool = pd.Series(duration_50_on_index(dur50, elev_index), index=elev_index)
        pools = {"rulecurve": float(rc_pool.loc[start]),
                 "duration50": float(d50_pool.loc[start]),
                 "observed": float(obs_pool.loc[start])
                 if np.isfinite(obs_pool.loc[start]) else np.nan}
        pools["highest"] = float(np.nanmax([pools["rulecurve"], pools["duration50"],
                                            pools["observed"]]))
        series_by_basis = {"rulecurve": rc_pool, "duration50": d50_pool,
                           "observed": obs_pool}
        # "highest" holds a constant at the highest of the three from the start
        series_by_basis["highest"] = pd.Series(pools["highest"], index=elev_index)
        events.append({"label": label, "note": note, "peak_date": centre,
                       "start": start, "index": index, "elev_index": elev_index,
                       "mos": mos_w, "cas": cas_w, "total": total,
                       "obs_peak": float(total.max()),
                       "local_share": float(cas_w.loc[total.idxmax()] /
                                            total.max()) if total.max() > 0 else np.nan,
                       "pools": pools, "pool_series": series_by_basis,
                       "missing": int(mos_w.isna().sum() + cas_w.isna().sum())})

    mapping_rows = []
    member = 0
    with HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION) as dst:
        for ev in events:
            for target_label, aep, target_peak in TARGETS:
                factor = target_peak / ev["obs_peak"]
                for basis in bases:
                    pool = ev["pools"].get(basis, np.nan)
                    basis_used = basis
                    if not np.isfinite(pool):
                        # Dec 1933 predates the observed pool record (starts Oct
                        # 1973). Substitute the rule curve and record that the
                        # member is not on the basis its label implies.
                        pool = ev["pools"]["rulecurve"]
                        basis_used = "%s->rulecurve (no observed pool)" % basis
                    member += 1
                    f_part = "C:%06d|%s" % (member, ENS_SUFFIX)
                    mos_v = to_dss_values(ev["mos"].values * factor)
                    cas_v = to_dss_values(ev["cas"].values * factor)
                    pool_series = ev["pool_series"][basis]
                    if basis_used != basis:
                        pool_series = ev["pool_series"]["rulecurve"]
                    elev_v = to_dss_values(pool_series.values)
                    rc_v = to_dss_values(rule_curve_on_index(ev["elev_index"]))
                    for parts, vals, units, dpart, wstart in [
                            (("", "MOSSYROCK", "FLOW-IN"), mos_v, "CFS", flow_d, flow_write),
                            (("", "CASTLE ROCK", "FLOW-LOCAL"), cas_v, "CFS", flow_d, flow_write),
                            (("", "MOS", "ELEV"), elev_v, "FEET", elev_d, elev_write),
                            (("", "MOS", "ELEV-RULECURVE"), rc_v, "FEET", elev_d, elev_write)]:
                        pathname = "/%s/%s/%s/%s/1HOUR/%s/" % (
                            parts[0], parts[1], parts[2], dpart, f_part)
                        dst.put_ts(build_container(pathname, vals, wstart, units,
                                                   "INST-VAL", 60))
                    synth_year = SYNTH_YEAR_BASE + member
                    synth_start = ev["start"] + pd.DateOffset(
                        years=synth_year - ev["start"].year)
                    mapping_rows.append({
                        "member": member, "ensemble_f_part": f_part,
                        "event": ev["label"], "event_note": ev["note"],
                        "target": target_label, "target_aep": aep,
                        "target_unreg_peak_cfs": target_peak,
                        "scale_factor": round(factor, 4),
                        "observed_unreg_peak_cfs": round(ev["obs_peak"], 1),
                        "scaled_unreg_peak_cfs": round(ev["obs_peak"] * factor, 1),
                        "local_share_at_peak": round(ev["local_share"], 4),
                        "pool_basis": basis,
                        "pool_basis_used": basis_used,
                        "start_pool_ft": round(pool, 2),
                        "synth_water_year": synth_year,
                        "real_start": synth_start,
                        "real_end": synth_start + pd.Timedelta(hours=n_hours - 1),
                        "source_start": ev["start"],
                        "source_peak_date": ev["peak_date"],
                        "ensemble_start": ens_start, "hours": n_hours,
                        "elev_ensemble_start": elev_ens_start,
                        "elev_hours": elev_hours,
                        "missing_input_hours": ev["missing"]})

    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(MAPPING_CSV, index=False)
    plot_events(events, mapping, PLOT_STEM)
    plot_pools(events, PLOT_STEM)

    print("\nSOURCE EVENTS")
    for ev in events:
        print("   %-8s peak %7.0f  local %4.1f%%  gaps %d  pools: rc %.1f / d50 %.1f / obs %s"
              % (ev["label"], ev["obs_peak"], 100 * ev["local_share"], ev["missing"],
                 ev["pools"]["rulecurve"], ev["pools"]["duration50"],
                 "%.1f" % ev["pools"]["observed"]
                 if np.isfinite(ev["pools"]["observed"]) else "n/a"))
    print("\nSCALE FACTORS")
    pivot = mapping[mapping["pool_basis"] == bases[0]].pivot(
        index="event", columns="target", values="scale_factor")
    print(pivot.to_string())
    big = mapping[mapping["scale_factor"] > 2.0]
    if len(big):
        print("\n   NOTE: %d members scale by more than 2x. The further a factor is"
              % len(big))
        print("   from 1, the more the synthetic depends on the assumption that")
        print("   hydrograph shape is preserved with magnitude.")
    print("\nSynthetic water years %d-%d, one per member, so reassembled blocks"
          % (mapping["synth_water_year"].min(), mapping["synth_water_year"].max()))
    print("never overlap. source_start holds the true event date.")
    print("\nMembers written : %d   records: %d (4 per member)"
          % (len(mapping), len(mapping) * 4))
    print("Mapping CSV     : %s" % MAPPING_CSV)
    print("Plots           : %s_events.png, %s_pools.png" % (PLOT_STEM, PLOT_STEM))


main()
