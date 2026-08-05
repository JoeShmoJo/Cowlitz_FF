#MOS_Special_Release_MinFloodPool.py
# -*- coding: utf-8 -*-
"""
MOS Special Flood Release screening -- minimum flood pool start assumption.

For each analyst-defined event (Events_Detailed.csv), the pool is re-started at
the seasonal rule curve (minimum flood pool) elevation on the event's first day.
The observed elevation rise during the event is carried forward AS STORAGE, not
as elevation, so the flatter storage-elevation relationship low in the pool does
not get mis-applied high in the pool:

    S_shift(d) = S(rulecurve(d0)) + [ S(obs(d)) - S(obs(d0)) ]
    elev_shift(d) = elev( S_shift(d) )

The shifted elevation and the hourly inflow are then looked up in the MOS inflow
ESRD special-curve table to get the prescribed release at each timestep.  A
nonzero prescribed release means the project would have been in Special Flood
Releases.

The prescribed release is then compared, as a rolling 24-hour volume, against
the release the project actually made -- Mayfield outflow less the Mayfield
local when the local is available, and Mayfield outflow alone when it is not.
Hourly negatives in that difference are timing noise and are left alone; only a
negative 24-HOUR volume is treated as a red flag on the data.  If the observed 24-hour volume already
equals or exceeds the prescribed volume, the observed regulation was at least as
aggressive as a rule-curve start would have required, despite starting lower in
the pool.  Where the prescribed volume is larger, the project would have had to
release more.  The 24-hour window absorbs timing mismatches between the computed
inflow and the shaped local.

Outputs: hourly CSV, event summary CSV, a DSS file of every plotted series, and
four plots.
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
from scipy.interpolate import RegularGridInterpolator
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
DSS_FILE = r"../data/ObsData_RegUnreg.dss"
ESRD_CSV = r"../data/MOS_Inflow_ESRD.csv"
STORAGE_CSV = r"../data/MOS_stor_rating_curve.csv"
EVENTS_CSV = r"../data/Events_Detailed.csv"
OUT_DIR = r"../output"
OUT_DSS = r"../output/MOS_Special_Release.dss"

PATH_RC = "//MOS/ELEV-RULECURVE/13Jul2012 - 19May2026/1Hour/CENWP-CALC/"
PATH_ELEV = "//MOS/ELEV-USGS/30Sep1973 - 16May2026/1Day/USGS/"
PATH_FLOW = "//MOS/FLOW-IN-CALC-CLEANED-VOLCOR/01Sep2008 - 01Apr2026/1Hour/CWMS/"
PATH_MAY_LOCAL = "//MAY/FLOW-LOCAL-SHAPED/28Oct2008 - 29Mar2026/1Hour/CWMS/"
PATH_MAY_OUT = ("/COWLITZ RIVER BELOW MAYFIELD DAM, WA/14238000/FLOW/"
                "01Jan1900 - 01Jan2100/1Hour/USGS/")

DSS_F_PART = "ESRD-SCREEN"      # F-part applied to every series written out
VOLUME_WINDOW_HOURS = 24        # window for the downstream volume comparison
VOLUME_CENTERED = True          # centered window absorbs timing error both ways
CLIP_NEGATIVE_RELEASE = False   # hourly negatives are timing noise; judge on the 24-hr volume
MIN_CONTROLLING_RELEASE = 5000.0  # project minimum release; below this the ESRD is not controlling
ELEV_GAP_FILL_DAYS = 5          # max internal gap in daily elevation to interpolate
MAX_POOL_ELEV = 778.5           # normal full pool at Mossyrock
CAP_AT_MAX_POOL = False         # True = hold the shifted pool at MAX_POOL_ELEV
STORAGE_EXTRAP_TO = 800.0       # extend the storage curve this high (flag, not clip)
WATER_YEAR_START_MONTH = 10

DSS_UNDEFINED = -3.4028234663852886e38
CFSHR_TO_ACREFT = 3600.0 / 43560.0

# ----------------------------------------------------------------------------


def read_dss_series(dss_file, pathname):
    """Read a DSS regular time series into a pandas Series (DSS end-of-period stamps)."""
    dss = HecDss.Open(dss_file)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        nodata = np.array(ts.nodata, dtype=bool)
        values[nodata] = np.nan
        values[values <= -900.0] = np.nan
        index = pd.DatetimeIndex([h.datetime() for h in ts.times])
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index()


def to_calendar_hourly(series):
    """Hour-ending stamps -> hour-beginning stamps (calendar attribution)."""
    out = series.copy()
    out.index = out.index - pd.Timedelta(hours=1)
    return out


def to_calendar_daily(series):
    """1DAY end-of-day stamps (midnight of the next day) -> calendar day stamps."""
    out = series.copy()
    out.index = (out.index - pd.Timedelta(days=1)).normalize()
    return out


def fill_short_gaps(series, max_days):
    """Linearly interpolate internal gaps in a daily series up to max_days long."""
    full = series.asfreq("D")
    return full.interpolate(method="time", limit=max_days, limit_area="inside")


def build_seasonal_rule_curve(rc_hourly):
    """Collapse the rule curve to a repeating (month, day) daily climatology."""
    daily = rc_hourly.resample("D").mean().dropna()
    return daily.groupby([daily.index.month, daily.index.day]).mean()


def rule_curve_on_dates(seasonal, dates):
    """Look up the seasonal rule curve for a DatetimeIndex (Feb 29 -> Feb 28)."""
    out = np.full(len(dates), np.nan)
    for i in range(len(dates)):
        month, day = int(dates.month[i]), int(dates.day[i])
        if (month, day) in seasonal.index:
            out[i] = seasonal.loc[(month, day)]
        elif month == 2 and day == 29:
            out[i] = seasonal.loc[(2, 28)]
    return pd.Series(out, index=dates)


def load_storage_curve(csv_file, extrap_to):
    """Elevation-storage rating curve, extended above the top of the table."""
    table = pd.read_csv(csv_file)
    elev = table.iloc[:, 0].values.astype(float)
    stor = table.iloc[:, 1].values.astype(float)
    order = np.argsort(elev)
    elev, stor = elev[order], stor[order]
    if extrap_to > elev[-1]:
        slope = (stor[-1] - stor[-2]) / (elev[-1] - elev[-2])
        extra_elev = np.arange(elev[-1] + 1.0, extrap_to + 0.5, 1.0)
        extra_stor = stor[-1] + slope * (extra_elev - elev[-1])
        elev = np.concatenate((elev, extra_elev))
        stor = np.concatenate((stor, extra_stor))
    return elev, stor


def elev_to_storage(elev_curve, stor_curve, elevation):
    return np.interp(np.asarray(elevation, dtype=float), elev_curve, stor_curve,
                     left=np.nan, right=np.nan)


def storage_to_elev(elev_curve, stor_curve, storage):
    return np.interp(np.asarray(storage, dtype=float), stor_curve, elev_curve,
                     left=np.nan, right=np.nan)


def load_esrd_table(csv_file):
    """Load the inflow ESRD special-curve table; returns elev grid, inflow grid, releases."""
    table = pd.read_csv(csv_file)
    elev_grid = table["Elevation_ft"].values.astype(float)
    inflow_grid = table.columns[1:].astype(float).values
    release_grid = table.iloc[:, 1:].values.astype(float)
    inflow_grid = np.concatenate(([0.0], inflow_grid))
    release_grid = np.column_stack((np.zeros(len(elev_grid)), release_grid))
    return elev_grid, inflow_grid, release_grid


def make_esrd_lookup(elev_grid, inflow_grid, release_grid):
    return RegularGridInterpolator((elev_grid, inflow_grid), release_grid,
                                   method="linear", bounds_error=False, fill_value=None)


def esrd_release(interp, elevation, inflow, elev_grid, inflow_grid):
    """Prescribed release for arrays of elevation and inflow, clipped to the table extent."""
    elev = np.clip(np.asarray(elevation, dtype=float), elev_grid.min(), elev_grid.max())
    flow = np.clip(np.asarray(inflow, dtype=float), inflow_grid.min(), inflow_grid.max())
    good = np.isfinite(elev) & np.isfinite(flow)
    out = np.full(elev.shape, np.nan)
    if good.any():
        out[good] = interp(np.column_stack((elev[good], flow[good])))
    return np.where(np.isfinite(out), np.maximum(out, 0.0), np.nan)


def read_event_definitions(csv_file):
    """Analyst-defined event windows; End is inclusive of the whole calendar day."""
    table = pd.read_csv(csv_file, encoding="utf-8-sig")
    table["Start"] = pd.to_datetime(table["Start"], format="%d-%b-%y")
    table["End"] = pd.to_datetime(table["End"], format="%d-%b-%y")
    events = []
    for i, row in table.iterrows():
        events.append((i + 1, row["Start"], row["End"] + pd.Timedelta(hours=23)))
    return events


def water_year(timestamp):
    return timestamp.year + (1 if timestamp.month >= WATER_YEAR_START_MONTH else 0)


def rolling_volume(series, hours, centered):
    """Rolling mean flow over the window, and the equivalent volume in acre-feet."""
    mean_flow = series.rolling(hours, center=centered, min_periods=1).mean()
    return mean_flow, mean_flow * hours * CFSHR_TO_ACREFT


def build_event_frame(event_id, start, end, flow_hourly, elev_daily, rc_seasonal,
                      may_out, may_local, elev_curve, stor_curve):
    """Hourly frame for one event with the storage-based minimum flood pool shift."""
    inflow = flow_hourly.loc[start:end]
    if len(inflow) == 0 or inflow.isna().all():
        return None
    days = pd.DatetimeIndex(inflow.index.normalize().unique())

    obs_elev = elev_daily.reindex(days)
    rc_elev = rule_curve_on_dates(rc_seasonal, days)
    day0 = days[0]

    obs_stor = pd.Series(elev_to_storage(elev_curve, stor_curve, obs_elev.values), index=days)
    start_stor = float(elev_to_storage(elev_curve, stor_curve, rc_elev.loc[day0]))
    shift_stor = start_stor + (obs_stor - obs_stor.loc[day0])
    shift_elev = pd.Series(storage_to_elev(elev_curve, stor_curve, shift_stor.values), index=days)
    if CAP_AT_MAX_POOL:
        shift_elev = shift_elev.clip(upper=MAX_POOL_ELEV)
        shift_stor = pd.Series(elev_to_storage(elev_curve, stor_curve, shift_elev.values),
                               index=days)

    frame = pd.DataFrame({"inflow_cfs": inflow.values}, index=inflow.index)
    frame["event"] = event_id
    frame["day"] = frame.index.normalize()
    frame["elev_obs_ft"] = frame["day"].map(obs_elev)
    frame["elev_rulecurve_ft"] = frame["day"].map(rc_elev)
    frame["elev_shifted_ft"] = frame["day"].map(shift_elev)
    frame["stor_obs_acft"] = frame["day"].map(obs_stor)
    frame["stor_shifted_acft"] = frame["day"].map(shift_stor)
    frame["elev_offset_day1_ft"] = rc_elev.loc[day0] - obs_elev.loc[day0]
    frame["may_out_cfs"] = may_out.reindex(frame.index)
    frame["may_local_cfs"] = may_local.reindex(frame.index)
    frame["local_available"] = frame["may_local_cfs"].notna()
    raw = frame["may_out_cfs"] - frame["may_local_cfs"].fillna(0.0)
    frame["release_obs_raw_cfs"] = raw
    frame["release_obs_cfs"] = raw.clip(lower=0.0) if CLIP_NEGATIVE_RELEASE else raw
    return frame


def add_release_and_volumes(hourly, interp, elev_grid, inflow_grid):
    """ESRD releases and the 24-hour volume comparison against downstream flow."""
    hourly["release_shifted_cfs"] = esrd_release(
        interp, hourly["elev_shifted_ft"].values, hourly["inflow_cfs"].values,
        elev_grid, inflow_grid)
    hourly["release_observed_cfs"] = esrd_release(
        interp, hourly["elev_obs_ft"].values, hourly["inflow_cfs"].values,
        elev_grid, inflow_grid)
    hourly["in_special_release"] = hourly["release_shifted_cfs"] >= MIN_CONTROLLING_RELEASE

    pieces = []
    for _, block in hourly.groupby("event", sort=False):
        block = block.copy()
        rel_mean, rel_vol = rolling_volume(block["release_shifted_cfs"],
                                           VOLUME_WINDOW_HOURS, VOLUME_CENTERED)
        obs_mean, obs_vol = rolling_volume(block["release_obs_cfs"],
                                           VOLUME_WINDOW_HOURS, VOLUME_CENTERED)
        block["release_24hr_cfs"] = rel_mean
        block["release_24hr_acft"] = rel_vol
        block["release_obs_24hr_cfs"] = obs_mean
        block["release_obs_24hr_acft"] = obs_vol
        block["release_deficit_24hr_acft"] = rel_vol - obs_vol
        active = block["in_special_release"].rolling(
            VOLUME_WINDOW_HOURS, center=VOLUME_CENTERED, min_periods=1).max().astype(bool)
        block["observed_meets_prescribed"] = active & (obs_vol >= rel_vol)
        block["prescribed_exceeds_observed"] = active & (rel_vol > obs_vol)
        block["obs_volume_negative"] = obs_vol < 0
        # calendar-day blocks, as a cross-check on the rolling window
        day_rel = block.groupby("day")["release_shifted_cfs"].transform("mean") * 24 * CFSHR_TO_ACREFT
        day_obs = block.groupby("day")["release_obs_cfs"].transform("mean") * 24 * CFSHR_TO_ACREFT
        day_active = block.groupby("day")["in_special_release"].transform("max").astype(bool)
        block["release_daily_acft"] = day_rel
        block["release_obs_daily_acft"] = day_obs
        block["prescribed_exceeds_observed_daily"] = day_active & (day_rel > day_obs)
        block["prescribed_exceeds_observed_clean"] = (
            block["prescribed_exceeds_observed"] & ~block["obs_volume_negative"])
        pieces.append(block)
    return pd.concat(pieces)


def summarize_events(hourly, event_lookup):
    """One row per event summarizing the special-release screening."""
    rows = []
    for event_id, block in hourly.groupby("event", sort=False):
        triggered = block["in_special_release"]
        short = block["prescribed_exceeds_observed"]
        over_pool = block["elev_shifted_ft"] > MAX_POOL_ELEV
        rows.append({
            "event": event_id,
            "start": event_lookup[event_id][0],
            "end": event_lookup[event_id][1],
            "water_year": water_year(block.index[0]),
            "hours": len(block),
            "peak_inflow_cfs": block["inflow_cfs"].max(),
            "start_elev_obs_ft": block["elev_obs_ft"].iloc[0],
            "start_elev_rulecurve_ft": block["elev_rulecurve_ft"].iloc[0],
            "elev_offset_day1_ft": block["elev_offset_day1_ft"].iloc[0],
            "storage_rise_acft": block["stor_obs_acft"].max() - block["stor_obs_acft"].iloc[0],
            "max_elev_shifted_ft": block["elev_shifted_ft"].max(),
            "hours_above_max_pool": int(over_pool.sum()),
            "special_release": bool(triggered.any()),
            "first_special_release": block.index[triggered.values].min() if triggered.any() else pd.NaT,
            "hours_in_special": int(triggered.sum()),
            "peak_release_shifted_cfs": block["release_shifted_cfs"].max(),
            "peak_release_observed_cfs": block["release_observed_cfs"].max(),
            "special_release_observed":
                bool((block["release_observed_cfs"] >= MIN_CONTROLLING_RELEASE).any()),
            "local_available_pct": 100.0 * block["local_available"].mean(),
            "local_pct_of_inflow": 100.0 * block["may_local_cfs"].fillna(0.0).sum()
                / max(block["inflow_cfs"].sum(), 1.0),
            "hours_obs_volume_negative": int(block["obs_volume_negative"].sum()),
            "min_obs_release_24hr_acft": block["release_obs_24hr_acft"].min(),
            "days_prescribed_exceeds_observed_daily":
                int(block.loc[block["prescribed_exceeds_observed_daily"], "day"].nunique()),
            "max_daily_deficit_acft":
                (block["release_daily_acft"] - block["release_obs_daily_acft"]).max(),
            "peak_obs_release_cfs": block["release_obs_cfs"].max(),
            "max_release_24hr_acft": block["release_24hr_acft"].max(),
            "obs_release_24hr_at_peak_acft":
                block["release_obs_24hr_acft"].loc[block["release_24hr_acft"].idxmax()],
            "max_obs_release_24hr_acft": block["release_obs_24hr_acft"].max(),
            "max_deficit_24hr_acft": block["release_deficit_24hr_acft"].max(),
            "prescribed_exceeds_observed": bool(short.any()),
            "hours_prescribed_exceeds_observed": int(short.sum()),
            "hours_prescribed_exceeds_observed_clean":
                int(block["prescribed_exceeds_observed_clean"].sum()),
            "observed_meets_prescribed_all_hours":
                bool(block.loc[triggered.values, "observed_meets_prescribed"].all())
                if triggered.any() else True,
        })
    summary = pd.DataFrame(rows)
    summary["is_annual_max"] = False
    for _, group in summary.groupby("water_year"):
        summary.loc[group["peak_inflow_cfs"].idxmax(), "is_annual_max"] = True
    return summary


def hourly_to_por_series(hourly, column, por_index):
    """Scatter one event-scoped column onto the full hourly period of record."""
    series = pd.Series(np.nan, index=por_index)
    series.loc[hourly.index] = hourly[column].values
    return series


def write_dss_output(hourly, out_dss, por_index):
    """Write every plotted series to DSS as 1Hour regular records."""
    records = [
        ("//MOS/FLOW-IN-EVENT//1Hour/%s/" % DSS_F_PART, "inflow_cfs", "CFS", "INST-VAL"),
        ("//MOS/ELEV-OBS-EVENT//1Hour/%s/" % DSS_F_PART, "elev_obs_ft", "FEET", "INST-VAL"),
        ("//MOS/ELEV-RULECURVE-SEASONAL//1Hour/%s/" % DSS_F_PART, "elev_rulecurve_ft", "FEET", "INST-VAL"),
        ("//MOS/ELEV-SHIFTED-MINFLDPOOL//1Hour/%s/" % DSS_F_PART, "elev_shifted_ft", "FEET", "INST-VAL"),
        ("//MOS/STOR-SHIFTED-MINFLDPOOL//1Hour/%s/" % DSS_F_PART, "stor_shifted_acft", "ACRE-FT", "INST-VAL"),
        ("//MOS/FLOW-ESRD-SPECIAL-SHIFTED//1Hour/%s/" % DSS_F_PART, "release_shifted_cfs", "CFS", "INST-VAL"),
        ("//MOS/FLOW-ESRD-SPECIAL-OBS//1Hour/%s/" % DSS_F_PART, "release_observed_cfs", "CFS", "INST-VAL"),
        ("//MOS/FLOW-ESRD-SPECIAL-24HR//1Hour/%s/" % DSS_F_PART, "release_24hr_cfs", "CFS", "INST-VAL"),
        ("//MAY/FLOW-OUT-EVENT//1Hour/%s/" % DSS_F_PART, "may_out_cfs", "CFS", "INST-VAL"),
        ("//MAY/FLOW-LOCAL-EVENT//1Hour/%s/" % DSS_F_PART, "may_local_cfs", "CFS", "INST-VAL"),
        ("//MOS/FLOW-OUT-IMPLIED//1Hour/%s/" % DSS_F_PART, "release_obs_cfs", "CFS", "INST-VAL"),
        ("//MOS/FLOW-OUT-IMPLIED-24HR//1Hour/%s/" % DSS_F_PART, "release_obs_24hr_cfs", "CFS", "INST-VAL"),
        ("//MOS/FLOW-ESRD-DEFICIT-24HR//1Hour/%s/" % DSS_F_PART, "release_deficit_24hr_acft", "ACRE-FT", "INST-VAL"),
    ]
    if os.path.exists(out_dss):
        os.remove(out_dss)
    # DSS stamps are end-of-period; shift the hour-beginning index back to that convention
    start_time = (por_index[0] + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()
    dss = HecDss.Open(out_dss, mode="rw")
    try:
        for pathname, column, units, dtype in records:
            values = hourly_to_por_series(hourly, column, por_index).values.astype(float)
            values = np.where(np.isfinite(values), values, DSS_UNDEFINED)
            tsc = TimeSeriesContainer(pathname, len(values), 1, values=values.tolist(),
                                      start_time=start_time, data_units=units,
                                      data_type=dtype)
            dss.put_ts(tsc, store_flag=2)
    finally:
        dss.close()
    return [r[0] for r in records]


def plot_event_summary(summary, out_file):
    """Peak prescribed release, downstream comparison, and duration by event."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    x = np.arange(len(summary))
    labels = ["%s%s" % (t.strftime("%d %b %y"), "*" if a else "")
              for t, a in zip(summary["start"], summary["is_annual_max"])]

    ax = axes[0]
    ax.bar(x, summary["peak_inflow_cfs"], color="0.78", label="Peak event inflow")
    ax.bar(x, summary["peak_release_shifted_cfs"], color="#c0392b", width=0.55,
           label="Peak prescribed release (start at rule curve)")
    ax.bar(x, summary["peak_release_observed_cfs"], color="#2c7fb8", width=0.28,
           label="Peak prescribed release (observed start)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title("MOS Special Flood Releases -- events started from minimum flood pool "
                 "(storage-based shift).  * = water year maximum")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    ax.bar(x - 0.19, summary["max_release_24hr_acft"], width=0.38, color="#c0392b",
           label="Peak 24-hr prescribed release volume")
    ax.bar(x + 0.19, summary["obs_release_24hr_at_peak_acft"], width=0.38,
           color="#4c9a2a", label="Concurrent 24-hr observed MOS release (MAY out - local)")
    ax.plot(x, summary["max_obs_release_24hr_acft"], ls="none", marker="_",
            markersize=13, color="0.35", label="Max 24-hr observed release (any hour)")
    for i in np.where(summary["prescribed_exceeds_observed"].values)[0]:
        ax.annotate("exceeds", xy=(x[i], summary["max_release_24hr_acft"].iloc[i]),
                    xytext=(0, 4), textcoords="offset points", ha="center",
                    fontsize=7, color="k", fontweight="bold")  # prescribed > observed
    ax.set_ylabel("Volume (ac-ft)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[2]
    ax.bar(x, summary["hours_in_special"], color="#c0392b",
           label="Hours in special releases")
    ax.bar(x, summary["hours_prescribed_exceeds_observed"], color="k", width=0.35,
           label="Hours prescribed 24-hr volume exceeds observed release")
    ax.set_ylabel("Hours")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=7)
    for i in range(len(summary)):
        if summary["special_release"].iloc[i]:
            ax.get_xticklabels()[i].set_color("#c0392b")
            ax.get_xticklabels()[i].set_fontweight("bold")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_event_panels(hourly, summary, out_file):
    """One sub-plot per event that enters Special Flood Releases."""
    hits = summary[summary["special_release"]].reset_index(drop=True)
    if len(hits) == 0:
        return
    ncol = 3
    nrow = int(np.ceil(len(hits) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 3.7 * nrow), squeeze=False)

    for k in range(nrow * ncol):
        ax = axes[k // ncol][k % ncol]
        if k >= len(hits):
            ax.axis("off")
            continue
        row = hits.loc[k]
        block = hourly[hourly["event"] == row["event"]].sort_index()
        times = block.index

        ax.plot(times, block["inflow_cfs"], color="0.45", lw=1.1)
        ax.fill_between(times, 0, block["release_shifted_cfs"], color="#c0392b", alpha=0.5)
        ax.plot(times, block["release_24hr_cfs"], color="#c0392b", lw=1.6)
        ax.plot(times, block["release_obs_24hr_cfs"], color="#4c9a2a", lw=1.6)
        ax.axhline(0.0, color="0.6", lw=0.7)
        top = ax.get_ylim()[1]
        exceed = block["prescribed_exceeds_observed"].values.astype(bool)
        if exceed.any():
            ax.fill_between(times, ax.get_ylim()[0], top, where=exceed,
                            color="k", alpha=0.10, step="mid")
        flag = block["obs_volume_negative"].values.astype(bool)
        if flag.any():
            ax.fill_between(times, ax.get_ylim()[0], top, where=flag,
                            color="#e67e22", alpha=0.16, step="mid")
        ax.set_ylabel("Flow (cfs)", fontsize=8)
        ax.set_xlabel("")
        ax.tick_params(labelsize=7)
        ax.set_title("%s -- WY%d%s\nday-1 shift %+.1f ft,  local %.0f%% of inflow"
                     % (row["start"].strftime("%d %b %Y"), row["water_year"],
                        "  (WY max)" if row["is_annual_max"] else "",
                        row["elev_offset_day1_ft"], row["local_pct_of_inflow"]),
                     fontsize=9)
        ax.grid(alpha=0.25)
        for label in ax.get_xticklabels():
            label.set_rotation(25)
            label.set_horizontalalignment("right")

        ax2 = ax.twinx()
        ax2.plot(times, block["elev_shifted_ft"], color="#c0392b", lw=1.3, ls="--")
        ax2.plot(times, block["elev_obs_ft"], color="#2c7fb8", lw=1.0, ls=":")
        ax2.plot(times, block["elev_rulecurve_ft"], color="#8e44ad", lw=1.0)
        ax2.axhline(MAX_POOL_ELEV, color="k", lw=0.8, ls="-.")
        ax2.set_ylabel("Elevation (ft)", fontsize=8)
        ax2.tick_params(labelsize=7)

    handles = [
        Line2D([], [], color="0.45", lw=1.2, label="Inflow"),
        Line2D([], [], color="#c0392b", lw=6, alpha=0.5, label="Prescribed special release (hourly)"),
        Line2D([], [], color="#c0392b", lw=1.6, label="Prescribed release, 24-hr mean"),
        Line2D([], [], color="#4c9a2a", lw=1.6, label="Observed MOS release (MAY out - local), 24-hr mean"),
        Line2D([], [], color="k", lw=6, alpha=0.10, label="Prescribed 24-hr volume exceeds observed"),
        Line2D([], [], color="#e67e22", lw=6, alpha=0.16, label="RED FLAG: 24-hr (MAY out - local) is negative"),
        Line2D([], [], color="#c0392b", lw=1.3, ls="--", label="Shifted elevation (rule curve start)"),
        Line2D([], [], color="#2c7fb8", lw=1.0, ls=":", label="Observed elevation"),
        Line2D([], [], color="#8e44ad", lw=1.0, label="Seasonal rule curve"),
        Line2D([], [], color="k", lw=0.8, ls="-.", label="Max pool %.1f ft" % MAX_POOL_ELEV),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8.5, frameon=False)
    fig.suptitle("Events entering Special Flood Releases from a minimum flood pool start", fontsize=12)
    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_esrd_trajectories(hourly, elev_grid, inflow_grid, release_grid, out_file):
    """Event paths through the ESRD table, over the special-release contours."""
    fig, ax = plt.subplots(figsize=(11, 8))
    qmax = min(90000.0, inflow_grid.max())
    col = inflow_grid <= qmax
    contours = ax.contour(inflow_grid[col], elev_grid, release_grid[:, col],
                          levels=[1, 2000, 5000, 10000, 20000, 30000, 45000, 60000],
                          colors="0.55", linewidths=0.9)
    ax.clabel(contours, fmt="%.0f", fontsize=7)
    zero = ax.contour(inflow_grid[col], elev_grid, release_grid[:, col],
                      levels=[0.5], colors="k", linewidths=2.0)
    ax.clabel(zero, fmt={0.5: "special release threshold"}, fontsize=8)

    for _, block in hourly.groupby("event", sort=False):
        triggered = (block["release_shifted_cfs"] > 0).any()
        ax.plot(block["inflow_cfs"], block["elev_shifted_ft"],
                color="#c0392b" if triggered else "#7fb3d5",
                lw=1.5 if triggered else 0.9, alpha=0.85 if triggered else 0.6)
        ax.plot(block["inflow_cfs"], block["elev_obs_ft"], color="0.6", lw=0.6, alpha=0.45)

    ax.axhspan(MAX_POOL_ELEV, 800, color="0.85", alpha=0.6, zorder=0)
    ax.axhline(MAX_POOL_ELEV, color="k", lw=1.0, ls="-.")
    ax.set_xlim(0, qmax)
    ax.set_ylim(690, 790)
    ax.set_xlabel("Inflow (cfs)")
    ax.set_ylabel("Pool elevation (ft)")
    ax.set_title("Event trajectories on the MOS inflow ESRD (storage-based shift)")
    handles = [
        Line2D([], [], color="#c0392b", lw=1.5, label="Shifted path -- enters special releases"),
        Line2D([], [], color="#7fb3d5", lw=1.0, label="Shifted path -- no special release"),
        Line2D([], [], color="0.6", lw=0.8, label="Observed path"),
        Line2D([], [], color="k", lw=2.0, label="Special release threshold"),
        Line2D([], [], color="k", lw=1.0, ls="-.", label="Max pool %.1f ft" % MAX_POOL_ELEV),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_volume_comparison(summary, out_file):
    """Max 24-hour prescribed release volume against the downstream volume."""
    fig, ax = plt.subplots(figsize=(8.5, 8))
    hit = summary["special_release"].values
    exceed = summary["prescribed_exceeds_observed"].values
    ax.scatter(summary.loc[~hit, "obs_release_24hr_at_peak_acft"],
               summary.loc[~hit, "max_release_24hr_acft"],
               s=45, facecolor="#7fb3d5", edgecolor="0.3", label="No special release")
    ax.scatter(summary.loc[hit & ~exceed, "obs_release_24hr_at_peak_acft"],
               summary.loc[hit & ~exceed, "max_release_24hr_acft"],
               s=70, facecolor="#c0392b", edgecolor="k", label="Observed release already meets it")
    ax.scatter(summary.loc[exceed, "obs_release_24hr_at_peak_acft"],
               summary.loc[exceed, "max_release_24hr_acft"],
               s=110, facecolor="#c0392b", edgecolor="k", marker="D",
               label="Prescribed exceeds observed release")
    for _, row in summary[hit].iterrows():
        ax.annotate(row["start"].strftime("%b %y"),
                    (row["obs_release_24hr_at_peak_acft"], row["max_release_24hr_acft"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    top = max(summary["obs_release_24hr_at_peak_acft"].max(),
              summary["max_release_24hr_acft"].max()) * 1.08
    ax.plot([0, top], [0, top], color="k", lw=1.0, ls="--")
    ax.text(top * 0.62, top * 0.66, "1:1", fontsize=9, rotation=45)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("Concurrent 24-hr observed MOS release, MAY out - local (ac-ft)")
    ax.set_ylabel("Max 24-hr prescribed special release volume (ac-ft)")
    ax.set_title("Peak 24-hr prescribed special release vs. the observed release actually made")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    rc_hourly = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_RC))
    elev_daily = fill_short_gaps(to_calendar_daily(read_dss_series(DSS_FILE, PATH_ELEV)),
                                 ELEV_GAP_FILL_DAYS)
    flow_hourly = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_FLOW))
    may_local = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_MAY_LOCAL))
    may_out = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_MAY_OUT)).dropna()

    rc_seasonal = build_seasonal_rule_curve(rc_hourly)
    elev_curve, stor_curve = load_storage_curve(STORAGE_CSV, STORAGE_EXTRAP_TO)
    elev_grid, inflow_grid, release_grid = load_esrd_table(ESRD_CSV)
    interp = make_esrd_lookup(elev_grid, inflow_grid, release_grid)

    events = read_event_definitions(EVENTS_CSV)
    event_lookup = {e[0]: (e[1], e[2]) for e in events}
    frames = []
    skipped = []
    for event_id, start, end in events:
        frame = build_event_frame(event_id, start, end, flow_hourly, elev_daily,
                                  rc_seasonal, may_out, may_local, elev_curve, stor_curve)
        if frame is None:
            skipped.append((event_id, start))
        else:
            frames.append(frame)
    hourly = pd.concat(frames)
    hourly = add_release_and_volumes(hourly, interp, elev_grid, inflow_grid)
    hourly.index.name = "datetime"

    summary = summarize_events(hourly, event_lookup)

    hourly.to_csv(os.path.join(OUT_DIR, "MOS_Special_Release_Hourly.csv"), float_format="%.2f")
    summary.to_csv(os.path.join(OUT_DIR, "MOS_Special_Release_Events.csv"),
                   index=False, float_format="%.2f")

    por_index = pd.date_range(flow_hourly.index[0], flow_hourly.index[-1], freq="h")
    written = write_dss_output(hourly, OUT_DSS, por_index)

    plot_event_summary(summary, os.path.join(OUT_DIR, "MOS_Special_Release_Summary.png"))
    plot_event_panels(hourly, summary, os.path.join(OUT_DIR, "MOS_Special_Release_WaterYears.png"))
    plot_esrd_trajectories(hourly, elev_grid, inflow_grid, release_grid,
                           os.path.join(OUT_DIR, "MOS_Special_Release_ESRD_Paths.png"))
    plot_volume_comparison(summary, os.path.join(OUT_DIR, "MOS_Special_Release_Volumes.png"))

    annual = summary[summary["is_annual_max"]]
    print("minimum controlling release: %.0f cfs" % MIN_CONTROLLING_RELEASE)
    print("events defined: %d   analyzed: %d   skipped (no inflow): %s"
          % (len(events), len(summary), skipped if skipped else "none"))
    print("entering special releases: %d of %d   (water year maxima only: %d of %d)"
          % (int(summary["special_release"].sum()), len(summary),
             int(annual["special_release"].sum()), len(annual)))
    print("entering special releases from the observed start: %d"
          % int(summary["special_release_observed"].sum()))
    print("prescribed 24-hr volume exceeds the observed MOS release: %d events"
          % int(summary["prescribed_exceeds_observed"].sum()))
    print("same test on fixed calendar-day blocks: %d events, %d days total"
          % (int((summary["days_prescribed_exceeds_observed_daily"] > 0).sum()),
             int(summary["days_prescribed_exceeds_observed_daily"].sum())))
    print("observed release already meets the prescribed volume at every triggered hour: %d events"
          % int(summary["observed_meets_prescribed_all_hours"].sum()))
    print("RED FLAG -- 24-hr volume of (MAY out - local) is negative: %d hours in %d events"
          % (int(summary["hours_obs_volume_negative"].sum()),
             int((summary["hours_obs_volume_negative"] > 0).sum())))
    print("   flagged hours with a non-negative 24-hr observed volume: %d of %d"
          % (int(summary["hours_prescribed_exceeds_observed_clean"].sum()),
             int(summary["hours_prescribed_exceeds_observed"].sum())))
    print("MAY local as a share of MOS inflow volume: median %.0f%%, max %.0f%% (event %s)"
          % (summary["local_pct_of_inflow"].median(), summary["local_pct_of_inflow"].max(),
             summary.loc[summary["local_pct_of_inflow"].idxmax(), "start"].strftime("%b %Y")))
    print("mean day-1 elevation shift: %+.1f ft   events above max pool: %d"
          % (summary["elev_offset_day1_ft"].mean(),
             int((summary["hours_above_max_pool"] > 0).sum())))
    print("MAY local availability: %.0f%% of event hours"
          % (100 * hourly["local_available"].mean()))
    print("DSS records written to %s:" % OUT_DSS)
    for pathname in written:
        print("   %s" % pathname)
    show = summary[["start", "water_year", "is_annual_max", "peak_inflow_cfs",
                    "elev_offset_day1_ft", "max_elev_shifted_ft", "hours_in_special",
                    "peak_release_shifted_cfs", "max_release_24hr_acft",
                    "obs_release_24hr_at_peak_acft", "max_deficit_24hr_acft",
                    "hours_prescribed_exceeds_observed",
                    "hours_prescribed_exceeds_observed_clean",
                    "hours_obs_volume_negative",
                    "days_prescribed_exceeds_observed_daily",
                    "max_daily_deficit_acft", "local_pct_of_inflow"]].copy()
    show["start"] = show["start"].dt.strftime("%Y-%m-%d")
    print(show.round(1).to_string(index=False))


main()
