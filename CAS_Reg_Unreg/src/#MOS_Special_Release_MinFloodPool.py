#MOS_Special_Release_MinFloodPool.py
# -*- coding: utf-8 -*-
"""
MOS Special Flood Release screening -- minimum flood pool start assumption.

For each contiguous event in the MOS cleaned/volume-corrected inflow record, the
observed daily pool elevation is shifted by a constant offset so the event's
first day starts at the seasonal rule curve (minimum flood pool) elevation.  The
shifted elevation and the hourly inflow are then looked up in the MOS inflow
ESRD special-curve table to get the prescribed release at each timestep.  A
nonzero prescribed release means the project would have been in Special Flood
Releases.

Outputs: hourly results CSV, event summary CSV, and three plots.
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

# ----------------------------------------------------------------------------
# PATHS AND SETTINGS
# ----------------------------------------------------------------------------
DSS_FILE = r"../data/ObsData_RegUnreg.dss"
ESRD_CSV = r"../data/MOS_Inflow_ESRD.csv"
OUT_DIR = r"../output"

PATH_RC = "//MOS/ELEV-RULECURVE/13Jul2012 - 19May2026/1Hour/CENWP-CALC/"
PATH_ELEV = "//MOS/ELEV-USGS/30Sep1973 - 16May2026/1Day/USGS/"
PATH_FLOW = "//MOS/FLOW-IN-CALC-CLEANED-VOLCOR/01Sep2008 - 01Apr2026/1Hour/CWMS/"

ELEV_GAP_FILL_DAYS = 5      # max internal gap in daily elevation to interpolate
SHIFT_UP_ONLY = False       # True = never shift the pool down below observed
MIN_EVENT_HOURS = 24        # ignore contiguous blocks shorter than this
TABLE_ELEV_MAX = 778.0      # top row of the ESRD table
MAX_POOL_ELEV = 778.5       # normal full pool at Mossyrock
CAP_AT_MAX_POOL = False     # True = hold the shifted pool at MAX_POOL_ELEV

# ----------------------------------------------------------------------------


def read_dss_series(dss_file, pathname):
    """Read a DSS regular time series into a pandas Series (DSS end-of-period stamps)."""
    from pydsstools.heclib.dss import HecDss
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
    """Linearly interpolate internal gaps up to max_days long."""
    full = series.asfreq("D")
    filled = full.interpolate(method="time", limit=max_days, limit_area="inside")
    return filled


def build_seasonal_rule_curve(rc_hourly):
    """Collapse the rule curve to a repeating (month, day) daily climatology."""
    daily = rc_hourly.resample("D").mean().dropna()
    seasonal = daily.groupby([daily.index.month, daily.index.day]).mean()
    return seasonal


def rule_curve_on_dates(seasonal, dates):
    """Look up the seasonal rule curve for a DatetimeIndex (Feb 29 -> Feb 28)."""
    months = dates.month.values
    days = dates.day.values
    out = np.full(len(dates), np.nan)
    for i in range(len(dates)):
        m, d = int(months[i]), int(days[i])
        if (m, d) in seasonal.index:
            out[i] = seasonal.loc[(m, d)]
        elif m == 2 and d == 29:
            out[i] = seasonal.loc[(2, 28)]
    return pd.Series(out, index=dates)


def load_esrd_table(csv_file):
    """Load the inflow ESRD special-curve table; returns elev grid, inflow grid, release grid."""
    table = pd.read_csv(csv_file)
    elev_grid = table["Elevation_ft"].values.astype(float)
    inflow_grid = table.columns[1:].astype(float).values
    release_grid = table.iloc[:, 1:].values.astype(float)
    # prepend a zero-inflow column so lookups below the first table column behave
    inflow_grid = np.concatenate(([0.0], inflow_grid))
    release_grid = np.column_stack((np.zeros(len(elev_grid)), release_grid))
    return elev_grid, inflow_grid, release_grid


def make_esrd_lookup(elev_grid, inflow_grid, release_grid):
    """Bilinear interpolator over the ESRD table."""
    return RegularGridInterpolator(
        (elev_grid, inflow_grid), release_grid,
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


def find_events(flow_hourly, min_hours):
    """Contiguous blocks of non-missing hourly inflow -> list of (start, end) timestamps."""
    valid = flow_hourly.notna().values.astype(int)
    edges = np.diff(np.concatenate(([0], valid, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    events = []
    for a, b in zip(starts, ends):
        if (b - a) >= min_hours:
            events.append((flow_hourly.index[a], flow_hourly.index[b - 1]))
    return events


def build_event_frame(flow_hourly, elev_daily, rc_seasonal, event_index, start, end):
    """Hourly frame for one event with observed, rule-curve and shifted elevations."""
    inflow = flow_hourly.loc[start:end]
    days = pd.DatetimeIndex(inflow.index.normalize().unique())
    obs = elev_daily.reindex(days)
    rc = rule_curve_on_dates(rc_seasonal, days)

    day0 = days[0]
    offset = rc.loc[day0] - obs.loc[day0]
    if SHIFT_UP_ONLY:
        offset = max(offset, 0.0)
    shifted = obs + offset
    if CAP_AT_MAX_POOL:
        shifted = shifted.clip(upper=MAX_POOL_ELEV)

    frame = pd.DataFrame({"inflow_cfs": inflow.values}, index=inflow.index)
    frame["event"] = event_index
    frame["day"] = frame.index.normalize()
    frame["elev_obs_ft"] = frame["day"].map(obs)
    frame["elev_rulecurve_ft"] = frame["day"].map(rc)
    frame["elev_shifted_ft"] = frame["day"].map(shifted)
    frame["offset_ft"] = offset
    return frame


def summarize_events(hourly):
    """One row per event summarizing the special-release screening."""
    rows = []
    for event_id, block in hourly.groupby("event"):
        triggered = block["release_shifted_cfs"] > 0
        triggered_obs = block["release_observed_cfs"] > 0
        first = block.index[triggered.values].min() if triggered.any() else pd.NaT
        over = block["elev_shifted_ft"] > MAX_POOL_ELEV
        first_over = block.index[over.values].min() if over.any() else pd.NaT
        rows.append({
            "event": event_id,
            "start": block.index[0],
            "end": block.index[-1],
            "hours": len(block),
            "peak_inflow_cfs": block["inflow_cfs"].max(),
            "start_elev_obs_ft": block["elev_obs_ft"].iloc[0],
            "start_elev_rulecurve_ft": block["elev_rulecurve_ft"].iloc[0],
            "offset_ft": block["offset_ft"].iloc[0],
            "max_elev_shifted_ft": block["elev_shifted_ft"].max(),
            "special_release": bool(triggered.any()),
            "first_special_release": first,
            "hours_in_special": int(triggered.sum()),
            "peak_release_shifted_cfs": block["release_shifted_cfs"].max(),
            "special_release_observed": bool(triggered_obs.any()),
            "hours_in_special_observed": int(triggered_obs.sum()),
            "peak_release_observed_cfs": block["release_observed_cfs"].max(),
            "hours_above_max_pool": int(over.sum()),
            "first_above_max_pool": first_over,
        })
    return pd.DataFrame(rows)


def plot_event_summary(summary, out_file):
    """Peak prescribed release and duration in special releases, by event."""
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)
    x = np.arange(len(summary))
    labels = [t.strftime("%b %Y") for t in summary["start"]]
    hit = summary["special_release"].values

    ax = axes[0]
    ax.bar(x, summary["peak_inflow_cfs"], color="0.75", label="Peak event inflow")
    ax.bar(x, summary["peak_release_shifted_cfs"], color="#c0392b", width=0.55,
           label="Peak prescribed release (start at rule curve)")
    ax.bar(x, summary["peak_release_observed_cfs"], color="#2c7fb8", width=0.28,
           label="Peak prescribed release (observed start)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title("MOS Special Flood Releases -- events started from minimum flood pool")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    ax.bar(x, summary["hours_in_special"], color="#c0392b", label="Hours in special releases (rule curve start)")
    ax.bar(x, summary["hours_in_special_observed"], color="#2c7fb8", width=0.45,
           label="Hours in special releases (observed start)")
    ax.set_ylabel("Hours")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    for i in range(len(summary)):
        if hit[i]:
            ax.get_xticklabels()[i].set_color("#c0392b")
            ax.get_xticklabels()[i].set_fontweight("bold")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_triggered_events(hourly, summary, out_file):
    """Per-event panels: inflow, prescribed release, and shifted vs observed elevation."""
    hits = summary[summary["special_release"]].reset_index(drop=True)
    if len(hits) == 0:
        return
    ncol = 3
    nrow = int(np.ceil(len(hits) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.4 * nrow), squeeze=False)

    for k in range(nrow * ncol):
        ax = axes[k // ncol][k % ncol]
        if k >= len(hits):
            ax.axis("off")
            continue
        event_id = hits.loc[k, "event"]
        block = hourly[hourly["event"] == event_id]
        hours = (block.index - block.index[0]).total_seconds() / 86400.0

        ax.plot(hours, block["inflow_cfs"], color="0.45", lw=1.2, label="Inflow")
        ax.fill_between(hours, 0, block["release_shifted_cfs"], color="#c0392b",
                        alpha=0.55, label="Prescribed release")
        ax.set_ylabel("Flow (cfs)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title("%s  (offset %+.1f ft)" %
                     (block.index[0].strftime("%d %b %Y"), hits.loc[k, "offset_ft"]),
                     fontsize=9)
        ax.grid(alpha=0.25)

        ax2 = ax.twinx()
        ax2.plot(hours, block["elev_shifted_ft"], color="#c0392b", lw=1.4, ls="--")
        ax2.plot(hours, block["elev_obs_ft"], color="#2c7fb8", lw=1.0, ls=":")
        ax2.plot(hours, block["elev_rulecurve_ft"], color="#2ca02c", lw=1.0)
        ax2.axhline(MAX_POOL_ELEV, color="k", lw=0.8, ls="-.")
        ax2.set_ylabel("Elevation (ft)", fontsize=8)
        ax2.tick_params(labelsize=7)
        if k % ncol == 0:
            ax.set_xlabel("Days into event", fontsize=8)

    handles = [
        Line2D([], [], color="0.45", lw=1.2, label="Inflow"),
        Line2D([], [], color="#c0392b", lw=6, alpha=0.55, label="Prescribed special release"),
        Line2D([], [], color="#c0392b", lw=1.4, ls="--", label="Shifted elevation (rule curve start)"),
        Line2D([], [], color="#2c7fb8", lw=1.0, ls=":", label="Observed elevation"),
        Line2D([], [], color="#2ca02c", lw=1.0, label="Seasonal rule curve"),
        Line2D([], [], color="k", lw=0.8, ls="-.", label="Max pool %.1f ft" % MAX_POOL_ELEV),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=9, frameon=False)
    fig.suptitle("Events entering Special Flood Releases from a minimum flood pool start", fontsize=12)
    fig.tight_layout(rect=[0, 0.045, 1, 0.97])
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_esrd_trajectories(hourly, elev_grid, inflow_grid, release_grid, out_file):
    """Event paths through the ESRD table, over the special-release contours."""
    fig, ax = plt.subplots(figsize=(11, 8))
    qmax = min(90000.0, inflow_grid.max())
    col = inflow_grid <= qmax
    levels = [1, 2000, 5000, 10000, 20000, 30000, 45000, 60000]
    cs = ax.contour(inflow_grid[col], elev_grid, release_grid[:, col],
                    levels=levels, colors="0.55", linewidths=0.9)
    ax.clabel(cs, fmt="%.0f", fontsize=7)
    zero = ax.contour(inflow_grid[col], elev_grid, release_grid[:, col],
                      levels=[0.5], colors="k", linewidths=2.0)
    ax.clabel(zero, fmt={0.5: "special release threshold"}, fontsize=8)

    for event_id, block in hourly.groupby("event"):
        triggered = (block["release_shifted_cfs"] > 0).any()
        ax.plot(block["inflow_cfs"], block["elev_shifted_ft"],
                color="#c0392b" if triggered else "#7fb3d5",
                lw=1.4 if triggered else 0.8, alpha=0.85 if triggered else 0.6)
        ax.plot(block["inflow_cfs"], block["elev_obs_ft"],
                color="0.6", lw=0.5, alpha=0.4)

    ax.axhspan(MAX_POOL_ELEV, 810, color="0.85", alpha=0.6, zorder=0)
    ax.axhline(MAX_POOL_ELEV, color="k", lw=1.0, ls="-.")
    ax.text(qmax * 0.985, MAX_POOL_ELEV + 1.5, "above max pool %.1f ft" % MAX_POOL_ELEV,
            ha="right", va="bottom", fontsize=8, color="0.35")
    ax.set_xlim(0, qmax)
    ax.set_ylim(680, 805)
    ax.set_xlabel("Inflow (cfs)")
    ax.set_ylabel("Pool elevation (ft)")
    ax.set_title("Event trajectories on the MOS inflow ESRD (special curves)")
    handles = [
        Line2D([], [], color="#c0392b", lw=1.4, label="Shifted event path -- enters special releases"),
        Line2D([], [], color="#7fb3d5", lw=1.0, label="Shifted event path -- no special release"),
        Line2D([], [], color="0.6", lw=0.8, label="Observed event path"),
        Line2D([], [], color="k", lw=2.0, label="Special release threshold"),
        Line2D([], [], color="k", lw=1.0, ls="-.", label="Max pool %.1f ft" % MAX_POOL_ELEV),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    rc_hourly = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_RC))
    elev_daily = to_calendar_daily(read_dss_series(DSS_FILE, PATH_ELEV))
    flow_hourly = to_calendar_hourly(read_dss_series(DSS_FILE, PATH_FLOW))

    elev_daily = fill_short_gaps(elev_daily, ELEV_GAP_FILL_DAYS)
    rc_seasonal = build_seasonal_rule_curve(rc_hourly)

    elev_grid, inflow_grid, release_grid = load_esrd_table(ESRD_CSV)
    interp = make_esrd_lookup(elev_grid, inflow_grid, release_grid)

    events = find_events(flow_hourly, MIN_EVENT_HOURS)
    frames = []
    for i, (start, end) in enumerate(events, start=1):
        frames.append(build_event_frame(flow_hourly, elev_daily, rc_seasonal, i, start, end))
    hourly = pd.concat(frames)

    hourly["release_shifted_cfs"] = esrd_release(
        interp, hourly["elev_shifted_ft"].values, hourly["inflow_cfs"].values,
        elev_grid, inflow_grid)
    hourly["release_observed_cfs"] = esrd_release(
        interp, hourly["elev_obs_ft"].values, hourly["inflow_cfs"].values,
        elev_grid, inflow_grid)
    hourly["in_special_release"] = hourly["release_shifted_cfs"] > 0
    hourly["elev_above_table"] = hourly["elev_shifted_ft"] > TABLE_ELEV_MAX
    hourly.index.name = "datetime"

    summary = summarize_events(hourly)

    hourly.to_csv(os.path.join(OUT_DIR, "MOS_Special_Release_Hourly.csv"),
                  float_format="%.2f")
    summary.to_csv(os.path.join(OUT_DIR, "MOS_Special_Release_Events.csv"),
                   index=False, float_format="%.2f")

    plot_event_summary(summary, os.path.join(OUT_DIR, "MOS_Special_Release_Summary.png"))
    plot_triggered_events(hourly, summary, os.path.join(OUT_DIR, "MOS_Special_Release_Events.png"))
    plot_esrd_trajectories(hourly, elev_grid, inflow_grid, release_grid,
                           os.path.join(OUT_DIR, "MOS_Special_Release_ESRD_Paths.png"))

    print("events screened: %d" % len(summary))
    print("events entering special releases (rule curve start): %d"
          % int(summary["special_release"].sum()))
    print("events entering special releases (observed start): %d"
          % int(summary["special_release_observed"].sum()))
    print("mean start offset: %+.1f ft   (negative offsets: %d)"
          % (summary["offset_ft"].mean(), int((summary["offset_ft"] < 0).sum())))
    print("events whose shifted pool exceeds max pool (%.1f ft): %d"
          % (MAX_POOL_ELEV, int((summary["hours_above_max_pool"] > 0).sum())))
    print("hours clipped at table top (%.0f ft): %d"
          % (TABLE_ELEV_MAX, int(hourly["elev_above_table"].sum())))
    show = summary[["start", "peak_inflow_cfs", "offset_ft", "max_elev_shifted_ft",
                    "hours_in_special", "peak_release_shifted_cfs",
                    "hours_above_max_pool"]].copy()
    show["start"] = show["start"].dt.strftime("%Y-%m-%d")
    print(show.round(1).to_string(index=False))


main()
