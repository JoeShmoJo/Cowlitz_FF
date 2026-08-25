#Coweeman_Event_Plotly.py
# -*- coding: utf-8 -*-
"""
Interactive (plotly) view of the Coweeman event analysis: Coweeman flow
against BOTH the regulated and unregulated Cowlitz at Castle Rock, over the
period where all three overlap, with every peak the analysis uses marked.

WHY THIS EXISTS
    #Coweeman_Proportion.py and #Coweeman_RegPeak_Timing.py both reduce this
    record to a table of ratios. The ratios are what feed the coincident
    frequency work, but they hide the thing you actually want to eyeball --
    whether the Coweeman is genuinely peaking alongside the Cowlitz, how far
    the regulated crest is displaced from the unregulated one, and whether a
    given event's numbers are trustworthy or an artifact of a gap in the
    Ecology record. A static PNG cannot support that; you have to zoom.

    Four quantities are marked on every event, and the last one is the whole
    point:

        unreg peak        the event definition -- what the events were
                          selected on, and what the frequency curve is
                          indexed by
        reg peak          the crest the tributary is actually being added
                          to, displaced later and lower by Mossyrock/Riffe
        Coweeman peak     the tributary's own crest, wherever it falls
        Coweeman @ reg    the Coweeman flow AT THE HOUR OF THE REGULATED
                          CREST -- the numerator of the LagScalingFactor,
                          and the only Coweeman number that belongs in a
                          sum with the regulated peak

    LagScalingFactor = (Coweeman @ reg) / (Coweeman peak). Seeing those two
    markers separate on screen is the argument for why the factor is not 1.0.

WHAT IS PLOTTED, AND AT WHAT RESOLUTION
    Overview  : hourly, whole overlap period, all events marked. Coweeman is
                resampled with max() (not mean) so a 15-minute crest is not
                averaged away before you see it.
    Per-event : native resolution -- Coweeman 15-minute, Cowlitz hourly --
                for the events above FACET_MIN_CFS. The flashiness of the
                Coweeman against the routed Cowlitz is the reason not to
                resample here.

    Both default to a log y-axis, with a linear toggle. Log is the default on
    purpose: the Coweeman runs ~20x smaller than the Cowlitz, so on a linear
    axis it is a flat line at the bottom of the frame. It is also the correct
    scale for reading a RATIO, which is what this whole analysis is about --
    a constant ratio is a constant vertical offset on a log axis.

EVENT SELECTION
    Identical parameters to #Coweeman_Proportion.py (N_EVENTS, separation,
    floor, window) so the events on this plot are the same 79 rows that are
    in coweeman_proportion.csv. Do not tune them here independently -- if
    they drift, the plot stops describing the table.

READ THE RECORD'S LIMITS BEFORE READING THE PLOT
    The overlap is Nov 2006 - Apr 2019, thirteen water years, from Ecology
    gage 26C075. The 80 markers are storm events (~6 per water year), NOT 80
    years of record. The largest concurrent event, Nov 2006 at ~155,000 cfs,
    is about a 1-in-59 on the unregulated curve, and only three events in the
    whole record exceed 100,000 cfs. Any tail ratio taken off this plot is an
    extrapolation, not a measurement.

    THE RATING CEILING IS THE BINDING LIMIT, and it is drawn on the plot.
    Gage 26C075's rating tops out near 3,400 cfs. Above that Ecology reports
    no discharge at all, only quality code 254 -- 305 readings over 10 days in
    WY2016-WY2019, shaded red here. Two of those blocks sit across the crest
    of an event in this set:

        Dec 2015   16.2 h censored, over the 2nd largest concurrent event
        Mar 2017   19.8 h censored, across the regulated crest exactly

    Inside a red band the Coweeman was HIGHER than anything the gage can
    report, so the trace flatlines at its last in-rating value and the event
    peak shown is a LOWER BOUND, not a measurement. Dropping those events
    biases the sample low, because they are precisely the events where the
    tributary was largest. Treat them as censored, not as data.

INPUTS
    ../data/coweeman/*_FM.txt        Ecology 26C075, 15-minute (cached by
                                     #Coweeman_Timing.py)
    ../output/ResSim_WCM_RC.dss      Flow-UNREG and Flow at CastleRock_NWS

OUTPUTS
    ../output/diagnostics/coweeman_events_overview.html
    ../output/diagnostics/coweeman_events_detail.html

REQUIRES
    plotly (pip install plotly). Writes self-contained HTML with plotly.js
    inlined, so the files open offline and can be handed to someone who does
    not have python.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import glob
import re
from datetime import datetime

import numpy as np
import pandas as pd
from pydsstools.heclib.dss import HecDss

import sys
REPO_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "Modules"))
# The Ecology parser lives in /Modules because four scripts once carried
# copy-pasted copies of it and all four shared the same bug: quality code 254
# ("Rating table exceeded, data will not be reported") was parsed as a
# discharge of 254 cfs. See Modules/ecology_io.py.
from ecology_io import (read_ecology_cache, resample_censor_aware,
                        MISSING_CODES, TRUSTED_CODES, CODE_MEANING)
from ecology_io import censored_spans

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    raise SystemExit(
        "plotly is not installed.  pip install plotly\n"
        "It is only needed to BUILD these files -- the HTML it writes is "
        "self-contained and opens without it.")

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
CACHE_DIR = r"../data/coweeman"
DSS_PATH = r"../output/ResSim_WCM_RC.dss"
UNREG_PATHNAME = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"
REG_PATHNAME = "//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/"

OUT_DIR = r"../output/diagnostics"
OVERVIEW_HTML = os.path.join(OUT_DIR, "coweeman_events_overview.html")
DETAIL_HTML = os.path.join(OUT_DIR, "coweeman_events_detail.html")

# Event selection -- HELD IDENTICAL to #Coweeman_Proportion.py on purpose.
N_EVENTS = 80
EVENT_MIN_SEPARATION_DAYS = 7
MIN_EVENT_CFS = 20000.0
EVENT_WINDOW_HOURS = 48          # Coweeman peak sought within +/- this
REG_SEARCH_WINDOW_HOURS = 72     # regulated crest sought within +/- this

# Per-event facets. 60000 reproduces the ">60k" tail bin of the analysis
# (9 events, a 3x3 grid). Lower it to widen the panel set.
FACET_MIN_CFS = 60000.0
FACET_PAD_HOURS = 72             # plotted window either side of the unreg peak
FACET_COLS = 3

C_UNREG = "#7aa9d0"
C_REG = "#1a4f8a"
C_COW = "#b7410e"
C_COW_AT_REG = "#7b2d8e"

def first_stamp(ts):
    """First timestamp of a DSS record, tolerant of both pydsstools builds.

    ts.times yields HecTime here and plain strings on the user's Windows
    build; DSS also writes midnight as 24:00 of the previous day, which
    pd.Timestamp refuses. See CLAUDE.md.
    """
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    roll = False
    if " 24:" in text:
        text, roll = text.replace(" 24:", " 00:"), True
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M"):
        try:
            stamp = pd.Timestamp(datetime.strptime(text, fmt))
            return stamp + pd.Timedelta(days=1) if roll else stamp
        except ValueError:
            continue
    return pd.Timestamp(text)


def read_dss_series(path, pathname, label):
    """Hourly series out of DSS, D part blanked so every block is assembled.

    The catalog returns one pathname per storage block; reading a block
    specific path returns only that block. See CLAUDE.md.
    """
    version = 6 if open(path, "rb").read(16)[12] == 6 else 7
    dss = HecDss.Open(path, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        values[np.array(ts.nodata, dtype=bool)] = np.nan
        values[values <= -900.0] = np.nan
        step = pd.Timedelta(seconds=int(ts.interval))
        index = pd.date_range(first_stamp(ts) - step, periods=len(values),
                              freq=step)
    finally:
        dss.close()
    series = pd.Series(values, index=index).dropna().sort_index()
    print("   %-25s: %d hours, %s to %s"
          % (label, len(series), series.index.min().date(),
             series.index.max().date()))
    return series


def find_events(series, n, min_separation_days, floor):
    """Top n independent peaks above floor, separation enforced both ways."""
    working = series.copy()
    events, span = [], pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span)
                    & (working.index <= when + span)] = np.nan
    return sorted(events)


def window_slice(series, center, hours):
    lo = center - pd.Timedelta(hours=hours)
    hi = center + pd.Timedelta(hours=hours)
    return series.loc[(series.index >= lo) & (series.index <= hi)]


def value_at(series, when):
    """Exact reading at a timestamp, or the nearest within 30 minutes.

    Deliberately NOT an interpolation. #Coweeman_RegPeak_Timing.py found that
    resampling papered over a MISSING-flagged reading on 2017-03-16 and
    invented a value that was never gaged. If nothing is within half an hour
    this returns NaN and the event drops out of the marker set, which is the
    honest outcome.
    """
    if not len(series):
        return np.nan
    idx = series.index.get_indexer([when], method="nearest")[0]
    if idx < 0:
        return np.nan
    stamp = series.index[idx]
    if abs(stamp - when) > pd.Timedelta(minutes=30):
        return np.nan
    return float(series.iloc[idx])


def code_at(qual, when):
    """Ecology quality code nearest a timestamp, or 0 if nothing is close."""
    if not len(qual):
        return 0
    idx = qual.index.get_indexer([when], method="nearest")[0]
    if idx < 0 or abs(qual.index[idx] - when) > pd.Timedelta(minutes=30):
        return 0
    return int(qual.iloc[idx])


def build_events(cow, unreg, reg, qual):
    """One row per event, with all four marked quantities located."""
    rows = []
    for when, unreg_peak in find_events(unreg, N_EVENTS,
                                        EVENT_MIN_SEPARATION_DAYS,
                                        MIN_EVENT_CFS):
        cow_win = window_slice(cow, when, EVENT_WINDOW_HOURS)
        if not len(cow_win):
            continue
        reg_win = window_slice(reg, when, REG_SEARCH_WINDOW_HOURS)
        if not len(reg_win):
            continue

        cow_peak_time = cow_win.idxmax()
        cow_peak = float(cow_win.loc[cow_peak_time])
        reg_peak_time = reg_win.idxmax()
        reg_peak = float(reg_win.loc[reg_peak_time])
        cow_at_reg = value_at(cow, reg_peak_time)

        rows.append({
            "unreg_peak_time": when,
            "unreg_peak_cfs": unreg_peak,
            "reg_peak_time": reg_peak_time,
            "reg_peak_cfs": reg_peak,
            "cow_peak_time": cow_peak_time,
            "cow_peak_cfs": cow_peak,
            "cow_at_reg_cfs": cow_at_reg,
            "reg_lag_hours": (reg_peak_time - when).total_seconds() / 3600.0,
            "cow_lag_hours": (cow_peak_time - when).total_seconds() / 3600.0,
            "lag_factor": (cow_at_reg / cow_peak) if cow_peak > 0 else np.nan,
            "ratio_peak": cow_peak / unreg_peak,
            "qual_at_cow_peak": code_at(qual, cow_peak_time),
            "qual_at_reg_peak": code_at(qual, reg_peak_time),
        })
    table = pd.DataFrame(rows)
    print("\n   %d events located (>= %.0f cfs unregulated, %d-day separation)"
          % (len(table), MIN_EVENT_CFS, EVENT_MIN_SEPARATION_DAYS))
    return table


def log_linear_buttons():
    """Toggle shared by both figures. Log is the default -- see docstring."""
    return [dict(
        type="buttons", direction="right", showactive=True,
        x=0.0, xanchor="left", y=1.12, yanchor="bottom", pad=dict(b=4),
        buttons=[
            dict(label="log y", method="relayout",
                 args=[{"yaxis.type": "log"}]),
            dict(label="linear y", method="relayout",
                 args=[{"yaxis.type": "linear"}]),
        ])]


def marker_trace(times, values, name, color, symbol, size=9, extra=None):
    text = extra if extra is not None else ["" for _ in times]
    return go.Scatter(
        x=times, y=values, mode="markers", name=name,
        marker=dict(color=color, symbol=symbol, size=size,
                    line=dict(color="white", width=1)),
        text=text,
        hovertemplate="<b>%s</b><br>%%{x|%%Y-%%m-%%d %%H:%%M}"
                      "<br>%%{y:,.0f} cfs%%{text}<extra></extra>" % name)


def build_overview(cow, unreg, reg, events, qual):
    """Whole overlap period, hourly, every event marked."""
    # max() not mean() -- a 15-minute Coweeman crest must survive the
    # resample, otherwise the markers float above their own trace.
    cow_h = cow.resample("1h").max().dropna().round(1)
    # Flow is stored to ~6 decimal places by ResSim and written into the HTML
    # as text. Rounding to whole cfs is lossless for any purpose this plot
    # serves and roughly halves the file.
    unreg = unreg.round(0)
    reg = reg.round(0)

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=unreg.index, y=unreg.values, name="Cowlitz unregulated",
        line=dict(color=C_UNREG, width=1),
        hovertemplate="unreg<br>%{x|%Y-%m-%d %H:%M}<br>%{y:,.0f} cfs<extra></extra>"))
    fig.add_trace(go.Scattergl(
        x=reg.index, y=reg.values, name="Cowlitz regulated",
        line=dict(color=C_REG, width=1),
        hovertemplate="reg<br>%{x|%Y-%m-%d %H:%M}<br>%{y:,.0f} cfs<extra></extra>"))
    fig.add_trace(go.Scattergl(
        x=cow_h.index, y=cow_h.values, name="Coweeman (26C075)",
        line=dict(color=C_COW, width=1),
        hovertemplate="Coweeman<br>%{x|%Y-%m-%d %H:%M}<br>%{y:,.0f} cfs<extra></extra>"))
    # Gaps in this trace are not low water -- they are readings Ecology
    # declined to report. The red bands below mark them.

    fig.add_trace(marker_trace(
        events["unreg_peak_time"], events["unreg_peak_cfs"],
        "unreg peak", C_UNREG, "triangle-down", 10))
    fig.add_trace(marker_trace(
        events["reg_peak_time"], events["reg_peak_cfs"],
        "reg peak", C_REG, "square", 8,
        extra=["<br>%+.0f h from unreg peak" % v for v in events["reg_lag_hours"]]))
    fig.add_trace(marker_trace(
        events["cow_peak_time"], events["cow_peak_cfs"],
        "Coweeman peak", C_COW, "circle", 8,
        extra=["<br>%+.0f h from unreg peak<br>%.3f of unreg peak"
               % (l, r) for l, r in zip(events["cow_lag_hours"],
                                        events["ratio_peak"])]))
    at_reg = events.dropna(subset=["cow_at_reg_cfs"])
    fig.add_trace(marker_trace(
        at_reg["reg_peak_time"], at_reg["cow_at_reg_cfs"],
        "Coweeman @ reg peak", C_COW_AT_REG, "x", 9,
        extra=["<br>LagScalingFactor %.3f" % v for v in at_reg["lag_factor"]]))

    for a, b in censored_spans(qual):
        fig.add_vrect(x0=a, x1=b, fillcolor="#c62828", opacity=0.16,
                      line_width=0, layer="below")

    fig.update_layout(
        title=dict(text="Coweeman vs Cowlitz at Castle Rock — %d events, "
                        "%s to %s<br><sub>13 water years of overlap. Markers "
                        "are storm events, not annual peaks. Log axis: a "
                        "constant ratio is a constant vertical offset.</sub>"
                   % (len(events), events["unreg_peak_time"].min().date(),
                      events["unreg_peak_time"].max().date())),
        xaxis=dict(title="", rangeslider=dict(visible=True), type="date"),
        yaxis=dict(title="Flow (cfs)", type="log"),
        hovermode="closest", template="plotly_white",
        legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1, xanchor="right"),
        updatemenus=log_linear_buttons(), height=680, margin=dict(t=130))
    return fig


def build_detail(cow, unreg, reg, events, qual):
    """Native-resolution panel per large event."""
    picked = events[events["unreg_peak_cfs"] >= FACET_MIN_CFS].copy()
    picked = picked.sort_values("unreg_peak_cfs", ascending=False)
    if not len(picked):
        raise SystemExit("No events at or above FACET_MIN_CFS=%.0f" % FACET_MIN_CFS)

    rows = int(np.ceil(len(picked) / FACET_COLS))
    titles = ["%s — unreg %s cfs" % (r["unreg_peak_time"].strftime("%d %b %Y"),
                                     format(int(r["unreg_peak_cfs"]), ","))
              for _, r in picked.iterrows()]
    fig = make_subplots(rows=rows, cols=FACET_COLS, subplot_titles=titles,
                        vertical_spacing=0.09, horizontal_spacing=0.06)

    for i, (_, ev) in enumerate(picked.iterrows()):
        r, c = i // FACET_COLS + 1, i % FACET_COLS + 1
        show = (i == 0)
        centre = ev["unreg_peak_time"]
        u = window_slice(unreg, centre, FACET_PAD_HOURS)
        g = window_slice(reg, centre, FACET_PAD_HOURS)
        w = window_slice(cow, centre, FACET_PAD_HOURS)

        fig.add_trace(go.Scatter(
            x=u.index, y=u.values, name="Cowlitz unregulated",
            legendgroup="unreg", showlegend=show,
            line=dict(color=C_UNREG, width=1.6),
            hovertemplate="unreg %{y:,.0f} cfs<br>%{x|%d %b %H:%M}<extra></extra>"),
            row=r, col=c)
        fig.add_trace(go.Scatter(
            x=g.index, y=g.values, name="Cowlitz regulated",
            legendgroup="reg", showlegend=show,
            line=dict(color=C_REG, width=1.6),
            hovertemplate="reg %{y:,.0f} cfs<br>%{x|%d %b %H:%M}<extra></extra>"),
            row=r, col=c)
        fig.add_trace(go.Scatter(
            x=w.index, y=w.values, name="Coweeman (26C075)",
            legendgroup="cow", showlegend=show,
            line=dict(color=C_COW, width=1.4),
            hovertemplate="Coweeman %{y:,.0f} cfs<br>%{x|%d %b %H:%M}<extra></extra>"),
            row=r, col=c)

        for key_t, key_v, nm, col, sym, grp in [
                ("unreg_peak_time", "unreg_peak_cfs", "unreg peak", C_UNREG,
                 "triangle-down", "m_unreg"),
                ("reg_peak_time", "reg_peak_cfs", "reg peak", C_REG,
                 "square", "m_reg"),
                ("cow_peak_time", "cow_peak_cfs", "Coweeman peak", C_COW,
                 "circle", "m_cow")]:
            fig.add_trace(go.Scatter(
                x=[ev[key_t]], y=[ev[key_v]], mode="markers", name=nm,
                legendgroup=grp, showlegend=show,
                marker=dict(color=col, symbol=sym, size=11,
                            line=dict(color="white", width=1.2)),
                hovertemplate="%s<br>%%{y:,.0f} cfs<br>%%{x|%%d %%b %%H:%%M}"
                              "<extra></extra>" % nm),
                row=r, col=c)

        if np.isfinite(ev["cow_at_reg_cfs"]):
            fig.add_trace(go.Scatter(
                x=[ev["reg_peak_time"]], y=[ev["cow_at_reg_cfs"]],
                mode="markers", name="Coweeman @ reg peak",
                legendgroup="m_atreg", showlegend=show,
                marker=dict(color=C_COW_AT_REG, symbol="x", size=12,
                            line=dict(color="white", width=1.2)),
                hovertemplate="Coweeman @ reg peak<br>%%{y:,.0f} cfs"
                              "<br>LagScalingFactor %.3f<extra></extra>"
                              % ev["lag_factor"]),
                row=r, col=c)
            # The vertical gap between this and the Coweeman peak IS the
            # LagScalingFactor -- draw it so the eye lands on it.
            fig.add_shape(
                type="line", row=r, col=c,
                x0=ev["reg_peak_time"], x1=ev["reg_peak_time"],
                y0=ev["cow_at_reg_cfs"], y1=ev["cow_peak_cfs"],
                line=dict(color=C_COW_AT_REG, width=1, dash="dot"))

        for a, b in censored_spans(qual, w.index.min(), w.index.max()):
            fig.add_vrect(x0=a, x1=b, fillcolor="#c62828", opacity=0.18,
                          line_width=0, layer="below", row=r, col=c)

        fig.update_yaxes(type="log", row=r, col=c)

    fig.update_layout(
        title=dict(text="Coweeman event detail — native resolution "
                        "(Coweeman 15-min, Cowlitz hourly)<br><sub>Events at "
                        "or above %s cfs unregulated. Dotted drop = Coweeman "
                        "peak down to its value at the regulated crest, i.e. "
                        "the LagScalingFactor.</sub>"
                   % format(int(FACET_MIN_CFS), ",")),
        template="plotly_white", hovermode="x unified",
        height=max(420, 330 * rows), margin=dict(t=140),
        legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1,
                    xanchor="right"),
        updatemenus=log_linear_buttons())
    fig.update_yaxes(title_text="Flow (cfs)", col=1)
    return fig


def report(events):
    print("\n   LagScalingFactor = Coweeman at the regulated crest / Coweeman peak")
    binned = pd.cut(events["unreg_peak_cfs"], [20000, 40000, 60000, 1e9],
                    labels=["20-40k", "40-60k", ">60k"])
    print(events.groupby(binned, observed=True)[
        ["lag_factor", "ratio_peak", "reg_lag_hours"]]
        .agg(["count", "median"]).round(3).to_string())
    tail = events[events["unreg_peak_cfs"] >= FACET_MIN_CFS]
    print("\n   quality of the Coweeman reading behind each tail-bin peak "
          "(>= %s cfs unregulated):" % format(int(FACET_MIN_CFS), ","))
    for _, ev in tail.sort_values("unreg_peak_cfs", ascending=False).iterrows():
        code = ev["qual_at_cow_peak"]
        print("     %s  unreg %9s   Coweeman peak %7s   code %-3d %s"
              % (ev["unreg_peak_time"].strftime("%d %b %Y"),
                 format(int(ev["unreg_peak_cfs"]), ","),
                 format(int(ev["cow_peak_cfs"]), ","), code,
                 CODE_MEANING.get(code, "?")))
    clean = tail[tail["qual_at_cow_peak"].isin(TRUSTED_CODES)]
    print("     -> %d of %d tail-bin peaks rest on an in-rating gaged reading."
          % (len(clean), len(tail)))

    missing = int(events["cow_at_reg_cfs"].isna().sum())
    if missing:
        print("\n   %d event(s) carry no 'Coweeman @ reg' marker: the gage was "
              "above its rating (or absent) at the hour of the regulated "
              "crest, so the flow is UNKNOWN AND HIGH. It is not a zero, and "
              "the LagScalingFactor medians above are biased LOW by their "
              "absence." % missing)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Reading...")
    cow, qual = read_ecology_cache(CACHE_DIR)
    unreg = read_dss_series(DSS_PATH, UNREG_PATHNAME, "Castle Rock unregulated")
    reg = read_dss_series(DSS_PATH, REG_PATHNAME, "Castle Rock regulated")

    lo = max(cow.index.min(), unreg.index.min(), reg.index.min())
    hi = min(cow.index.max(), unreg.index.max(), reg.index.max())
    print("\n   overlap: %s to %s  (%.1f years)"
          % (lo.date(), hi.date(), (hi - lo).days / 365.25))
    cow = cow.loc[lo:hi]
    qual = qual.loc[lo:hi]
    unreg = unreg.loc[lo:hi]
    reg = reg.loc[lo:hi]

    events = build_events(cow, unreg, reg, qual)
    report(events)

    for fig, path, label in [
            (build_overview(cow, unreg, reg, events, qual), OVERVIEW_HTML, "overview"),
            (build_detail(cow, unreg, reg, events, qual), DETAIL_HTML, "detail")]:
        fig.write_html(path, include_plotlyjs=True, full_html=True)
        print("   wrote %s  (%s, %.1f MB)"
              % (path, label, os.path.getsize(path) / 1e6))


if __name__ == "__main__":
    main()
