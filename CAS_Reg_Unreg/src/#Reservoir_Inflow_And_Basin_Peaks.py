#Reservoir_Inflow_And_Basin_Peaks.py
# -*- coding: utf-8 -*-
"""
Reservoir Inflow Analysis Script
"""
#Reservoir_Inflow
from pydsstools.heclib.dss import HecDss

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# NOTE FOR REVIEWERS: this script WRITES records into ../../CAS_Unreg_FF/data/obsData.dss
# (the archival observed-data store used by the rest of the analysis).
# Re-running regenerates those records; do not run casually.

from pydsstools.core import TimeSeriesContainer
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import webbrowser
import os


def read_from_dss(dss_file, pathnames, startDate, endDate):
    DataDict = {}
    with HecDss.Open(dss_file) as fid:
        for pathname in pathnames:
            try:
                ts = fid.read_ts(pathname, window=(startDate, endDate))
                times = pd.to_datetime(ts.pytimes)
                values = ts.values
                df = pd.DataFrame({'value': values}, index=times)
                df[df <= -900] = np.nan
                DataDict[pathname] = df
                print(f"Read: {pathname}")
            except Exception as e:
                print(f"Failed to read {pathname}: {e}")
    return DataDict


def recalculate_inflow_variants(DataDict, reservoir,
                                stor_path_override=None,
                                outflow_path_override=None):
    """
    Calculate one inflow variant (raw, no smoothing).
    Stored under two keys:
      - FLOW-IN-CALC-RAW       -> written to DSS as full continuous record
      - FLOW-IN-CALC-RAW-PEAKS -> written to DSS only around peak windows
    Both contain identical data; the DSS write functions handle the difference.

    stor_path_override / outflow_path_override: if provided, use these
    DataDict keys instead of the default CWMS paths (e.g. for MAY cleaned paths).
    """
    ACFT_PER_CFS_HOUR = 3600 / 43560

    stor_path    = stor_path_override    or f'//{reservoir}/STOR//1HOUR/CWMS/'
    outflow_path = outflow_path_override or f'//{reservoir}/FLOW-OUT//1HOUR/CWMS/'

    stor_df    = DataDict.get(stor_path)
    outflow_df = DataDict.get(outflow_path)

    if stor_df is None or outflow_df is None:
        print(f"[WARNING] Missing data for {reservoir} "
              f"(stor={stor_path}, outflow={outflow_path}); skipping.")
        return DataDict

    stor    = stor_df['value'].copy().replace(-902, np.nan)
    outflow = outflow_df['value'].copy().replace(-902, np.nan)

    combined = pd.DataFrame({'stor': stor, 'outflow': outflow})
    combined['dS_cfs'] = combined['stor'].diff() / ACFT_PER_CFS_HOUR
    combined['inflow'] = combined['outflow'] + combined['dS_cfs']
    combined = combined.dropna(subset=['stor', 'outflow', 'dS_cfs'])
    raw_series = combined['inflow'].rename('value')

    raw_path       = f'//{reservoir}/FLOW-IN-CALC-RAW//1HOUR/CWMS/'
    raw_peaks_path = f'//{reservoir}/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/'

    DataDict[raw_path]       = raw_series.to_frame()
    DataDict[raw_peaks_path] = raw_series.to_frame()

    print(f"[INFO] Calculated {raw_path}: {len(raw_series)} values")
    print(f"[INFO] Calculated {raw_peaks_path}: {len(raw_series)} values")

    return DataDict


def get_water_year_peaks(series):
    """Return annual peak (date, value) per water year (Oct 1 - Sep 30)."""
    peaks = []
    series = series.dropna()
    wy = series.index.to_series().apply(
        lambda t: t.year if t.month >= 10 else t.year - 1
    )
    for year, group in series.groupby(wy):
        if group.empty:
            continue
        peak_idx = group.idxmax()
        peaks.append({'water_year': year, 'date': peak_idx, 'value': group[peak_idx]})
    return pd.DataFrame(peaks)


def compute_combined_peaks(ref_series, window_hours=36):
    """
    For each hourly timestep, sum the rolling max of each stream over a
    window_hours window. Find annual water year peaks of the combined series.
    """
    all_series = pd.DataFrame(ref_series).resample('h').mean()

    rolling_max = all_series.rolling(
        window=window_hours,
        center=True,
        min_periods=1
    ).max()

    combined_series = rolling_max.sum(axis=1).dropna()

    wy = combined_series.index.to_series().apply(
        lambda t: t.year if t.month >= 10 else t.year - 1
    )
    peaks = []
    for year, group in combined_series.groupby(wy):
        if group.empty:
            continue
        peak_idx = group.idxmax()
        peaks.append({
            'water_year':    year,
            'combined_date': peak_idx,
            'combined_flow': round(group[peak_idx], 0),
        })

    peaks_df = pd.DataFrame(peaks)
    peaks_df['combined_rank'] = peaks_df['combined_flow'].rank(
        ascending=False, method='min').astype(int)

    print(f"Combined peaks computed for {len(peaks_df)} water years.")
    return peaks_df


def write_full_record_to_dss(dss_file, DataDict, reservoir, pathname_suffix):
    """Write a full continuous record to DSS, filling gaps with sentinel -901."""
    pathname = f'//{reservoir}/{pathname_suffix}//1HOUR/CWMS/'
    df = DataDict.get(pathname)
    if df is None or df.empty:
        print(f"[WARNING] {pathname} not found; skipping full record write.")
        return

    series = df.iloc[:, 0].copy()
    full_index = pd.date_range(start=series.index.min(),
                               end=series.index.max(), freq='h')
    series = series.reindex(full_index).fillna(-901)

    tsc = TimeSeriesContainer()
    tsc.pathname      = pathname
    tsc.startDateTime = str(series.index[0])
    tsc.numberValues  = len(series)
    tsc.values        = series.to_numpy()
    tsc.interval      = 1
    tsc.units         = 'CFS'
    tsc.type          = 'INST-VAL'

    with HecDss.Open(dss_file, version=6) as fid:
        fid.put_ts(tsc)

    n_valid = (series > -900).sum()
    print(f"[INFO] Wrote full record: {pathname} ({n_valid} valid + "
          f"{len(series) - n_valid} sentinel values, "
          f"{series.index.min().date()} to {series.index.max().date()})")


def write_peak_windows_to_dss(dss_file, DataDict, reservoir, ref_paths,
                               combined_peaks, window_days=15):
    """
    Write FLOW-IN (original) and FLOW-IN-CALC-RAW-PEAKS to DSS only for
    periods around individual reference stream peaks and combined peaks.
    """
    window = pd.Timedelta(days=window_days)

    peak_dates = set()

    for label, path in ref_paths.items():
        df = DataDict.get(path)
        if df is None or df.empty:
            continue
        series = df.iloc[:, 0].replace(-902, np.nan).dropna()
        peaks_df = get_water_year_peaks(series)
        for _, row in peaks_df.iterrows():
            peak_dates.add(row['date'])

    for _, row in combined_peaks.iterrows():
        peak_dates.add(row['combined_date'])

    intervals = sorted((d - window, d + window) for d in peak_dates)
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    print(f"\n[INFO] Writing peak windows to DSS: {len(merged)} merged intervals "
          f"for {reservoir}")

    variant_paths = [
        f'//{reservoir}/FLOW-IN//1HOUR/CWMS/',
        f'//{reservoir}/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/',
    ]

    for pathname in variant_paths:
        df = DataDict.get(pathname)
        if df is None or df.empty:
            print(f"[WARNING] {pathname} not found; skipping.")
            continue

        series = df.iloc[:, 0].copy()

        for t_start, t_end in merged:
            windowed = series[(series.index >= t_start) & (series.index <= t_end)].dropna()
            if windowed.empty:
                continue

            tsc = TimeSeriesContainer()
            tsc.pathname      = pathname
            tsc.startDateTime = str(windowed.index[0])
            tsc.numberValues  = len(windowed)
            tsc.values        = windowed.to_numpy()
            tsc.interval      = 1
            tsc.units         = 'CFS'
            tsc.type          = 'INST-VAL'

            with HecDss.Open(dss_file, version=6) as fid:
                fid.put_ts(tsc)

            print(f"  Wrote {len(windowed)} values: {pathname} "
                  f"[{t_start.date()} to {t_end.date()}]")


def plot_peak_windows(DataDict, reservoir, ref_paths, combined_peaks,
                      window_days=15):
    """
    Single plot showing all water year peak windows for reference streams
    and reservoir inflow variants. Peaks marked with asterisk.
    Prints and saves a table of peak dates including combined peak info.
    """
    ref_series = {}
    for label, path in ref_paths.items():
        df = DataDict.get(path)
        if df is not None and not df.empty:
            ref_series[label] = df.iloc[:, 0].replace(-902, np.nan).dropna()

    if not ref_series:
        print("[WARNING] No reference stream data found.")
        return

    res_series = {
        f'{reservoir} Original':         DataDict.get(f'//{reservoir}/FLOW-IN//1HOUR/CWMS/'),
        f'{reservoir} Calc Raw':         DataDict.get(f'//{reservoir}/FLOW-IN-CALC-RAW//1HOUR/CWMS/'),
        f'{reservoir} Calc Raw (Peaks)': DataDict.get(f'//{reservoir}/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/'),
    }
    res_clean = {}
    for label, df in res_series.items():
        if df is not None and not df.empty:
            res_clean[label] = df.iloc[:, 0].replace(-902, np.nan)

    all_peaks = {}
    for label, series in ref_series.items():
        all_peaks[label] = get_water_year_peaks(series)
        print(f"{label}: {len(all_peaks[label])} peaks found")

    window_events = []
    for label, peaks_df in all_peaks.items():
        for _, peak_row in peaks_df.iterrows():
            window_events.append({
                'label':      label,
                'water_year': peak_row['water_year'],
                'date':       peak_row['date'],
                'value':      peak_row['value'],
            })

    combined_date_map = dict(zip(combined_peaks['water_year'], combined_peaks['combined_date']))
    combined_rank_map = dict(zip(combined_peaks['water_year'], combined_peaks['combined_rank']))
    combined_flow_map = dict(zip(combined_peaks['water_year'], combined_peaks['combined_flow']))

    peaks_table = pd.DataFrame(window_events).copy()
    peaks_table['date_str'] = peaks_table['date'].dt.strftime('%Y-%m-%d %H:%M')
    peaks_table['value'] = peaks_table['value'].round(0).astype(int)

    peaks_table['percentile'] = peaks_table.apply(
        lambda row: round(
            (sum(1 for e in window_events
                 if e['label'] == row['label'] and e['value'] <= row['value'])
             / sum(1 for e in window_events if e['label'] == row['label'])) * 100, 1
        ), axis=1
    )

    def days_from_combined(row):
        combined_date = combined_date_map.get(row['water_year'])
        if combined_date is None:
            return np.nan
        return round((row['date'] - combined_date).total_seconds() / 86400, 2)

    peaks_table['days_from_combined'] = peaks_table.apply(days_from_combined, axis=1)
    peaks_table['combined_date_str'] = pd.to_datetime(
        peaks_table['water_year'].map(lambda wy: combined_date_map.get(wy, pd.NaT))
    ).dt.strftime('%Y-%m-%d %H:%M')
    peaks_table['combined_flow'] = peaks_table['water_year'].map(
        lambda wy: combined_flow_map.get(wy, np.nan))
    peaks_table['combined_rank'] = peaks_table['water_year'].map(
        lambda wy: combined_rank_map.get(wy, np.nan))

    peaks_table_out = peaks_table[[
        'label', 'water_year', 'date_str', 'value', 'percentile',
        'days_from_combined', 'combined_date_str', 'combined_flow', 'combined_rank'
    ]].copy()
    peaks_table_out = peaks_table_out.sort_values(
        ['water_year', 'label']).reset_index(drop=True)
    peaks_table_out.columns = [
        'Stream', 'Water Year', 'Peak Date', 'Peak Flow (cfs)', 'Percentile',
        'Days from Combined Peak', 'Combined Peak Date',
        'Combined Flow (cfs)', 'Combined Rank'
    ]
    peaks_table_out.to_csv('../output/peak_dates.csv', index=False)
    print("\nPeak dates:")
    print(peaks_table_out.to_string(index=False))
    print(f"Saved: peak_dates.csv\n")

    ref_color_list = ['royalblue', 'darkorange', 'forestgreen', 'crimson',
                      'darkorchid', 'saddlebrown']
    ref_colors = {label: ref_color_list[i % len(ref_color_list)]
                  for i, label in enumerate(ref_paths.keys())}
    res_colors = {
        f'{reservoir} Original':         'black',
        f'{reservoir} Calc Raw':         'red',
        f'{reservoir} Calc Raw (Peaks)': 'darkred',
    }

    window = pd.Timedelta(days=window_days)
    fig = go.Figure()
    dummy_date = pd.Timestamp('1900-01-01')

    # Dummy traces for legend
    for label, color in ref_colors.items():
        fig.add_trace(go.Scatter(
            x=[dummy_date], y=[0], name=label,
            line=dict(color=color, width=1.5),
            legendgroup=label, showlegend=True, visible='legendonly',
        ))
        fig.add_trace(go.Scatter(
            x=[dummy_date], y=[0], name=f'{label} Peak',
            mode='markers',
            marker=dict(symbol='asterisk', size=14, color=color, line_width=2),
            legendgroup=f'{label}_peak', showlegend=True, visible='legendonly',
        ))

    for label, color in res_colors.items():
        if label not in res_clean:
            continue
        fig.add_trace(go.Scatter(
            x=[dummy_date], y=[0], name=label,
            line=dict(color=color, width=1,
                      dash='solid' if 'Original' in label else 'dot'),
            legendgroup=label, showlegend=True, visible='legendonly', opacity=0.8,
        ))

    # Plot windows
    for event in window_events:
        peak_date  = event['date']
        peak_label = event['label']
        t_start    = peak_date - window
        t_end      = peak_date + window

        for label, series in ref_series.items():
            windowed = series[(series.index >= t_start) & (series.index <= t_end)]
            if windowed.empty:
                continue
            fig.add_trace(go.Scatter(
                x=windowed.index, y=windowed.values, name=label,
                line=dict(color=ref_colors[label], width=1.5),
                legendgroup=label, showlegend=False, opacity=0.9,
            ))
            if label == peak_label and peak_date in series.index:
                fig.add_trace(go.Scatter(
                    x=[peak_date], y=[series[peak_date]],
                    mode='markers',
                    marker=dict(symbol='asterisk', size=14,
                                color=ref_colors[label], line_width=2),
                    name=f'{label} Peak',
                    legendgroup=f'{label}_peak', showlegend=False,
                ))

        for label, series in res_clean.items():
            windowed = series[(series.index >= t_start) & (series.index <= t_end)]
            if windowed.empty:
                continue
            fig.add_trace(go.Scatter(
                x=windowed.index, y=windowed.values, name=label,
                line=dict(color=res_colors[label], width=1,
                          dash='solid' if 'Original' in label else 'dot'),
                legendgroup=label, showlegend=False, opacity=0.8,
            ))

        fig.add_vline(
            x=peak_date,
            line=dict(color=ref_colors[peak_label], width=0.5, dash='dash'),
        )

    fig.update_layout(
        title=f'{reservoir} Inflow vs Unregulated Streams — ±{window_days} Day Windows Around WY Peaks',
        xaxis_title='Date',
        yaxis_title='Flow (cfs)',
        height=600,
        hovermode='x unified',
        legend=dict(groupclick='toggleitem'),
    )

    output_path = os.path.abspath(f'../diagnostics/{reservoir}_peak_windows.html')
    fig.write_html(output_path)
    print(f"Saved: {output_path}")
    webbrowser.open(f'file:///{output_path}')


def plot_inflow_interactive(DataDict, reservoirs):
    """
    Interactive Plotly plot of original vs raw calculated inflow.
    One subplot per reservoir, shared x-axis, all series toggleable.
    """
    fig = make_subplots(
        rows=len(reservoirs), cols=1,
        shared_xaxes=True,
        subplot_titles=[f'{r} — Inflow Comparison' for r in reservoirs],
        vertical_spacing=0.08
    )

    for i, reservoir in enumerate(reservoirs, start=1):
        series_map = {
            'Original (CWMS)':  (DataDict.get(f'//{reservoir}/FLOW-IN//1HOUR/CWMS/'),                 'steelblue'),
            'Calc Raw':         (DataDict.get(f'//{reservoir}/FLOW-IN-CALC-RAW//1HOUR/CWMS/'),        'red'),
            'Calc Raw (Peaks)': (DataDict.get(f'//{reservoir}/FLOW-IN-CALC-RAW-PEAKS//1HOUR/CWMS/'), 'darkred'),
        }

        for label, (df, color) in series_map.items():
            if df is None or df.empty:
                continue
            series = df.iloc[:, 0].replace(-902, np.nan)
            fig.add_trace(go.Scatter(
                x=series.index,
                y=series.values,
                name=f'{reservoir} {label}',
                line=dict(color=color, width=1),
                legendgroup=f'{reservoir}_{label}',
                opacity=0.85,
            ), row=i, col=1)

        fig.update_yaxes(title_text='Flow (cfs)', row=i, col=1)

    fig.update_layout(
        height=400 * len(reservoirs),
        title_text='Reservoir Inflow Comparison',
        title_font_size=15,
        hovermode='x unified',
        legend=dict(groupclick='toggleitem'),
    )
    fig.update_xaxes(title_text='Date', row=len(reservoirs), col=1)

    output_path = os.path.abspath('../diagnostics/inflow_comparison.html')
    fig.write_html(output_path)
    print(f"Saved: {output_path}")
    webbrowser.open(f'file:///{output_path}')


# =============================================================================
# MAIN
# =============================================================================

ref_paths = {
    'Cowlitz at Packwood': '/COWLITZ RIVER AT PACKWOOD, WA/14226500/FLOW//1HOUR/USGS/',
    'Tilton at Cinebar':   '/TILTON RIVER ABOVE BEAR CANYON CREEK NEAR CINEBAR, WA/14236200/FLOW//1HOUR/USGS/',
    'Toutle at Tower Rd':  '/TOUTLE RIVER AT TOWER ROAD NEAR SILVER LAKE, WA/14242580/FLOW//1HOUR/USGS/',
}

local_flow_paths = {
    'Cowlitz at Castle Rock':     '/COWLITZ RIVER AT CASTLE ROCK, WA/14243000/FLOW/01Oct2008 - 11May2026/1HOUR/USGS/',
    'Cowlitz below Mayfield Dam': '/COWLITZ RIVER BELOW MAYFIELD DAM, WA/14238000/FLOW/01Oct2008 - 11May2026/1HOUR/USGS/',
}

# MAY-specific cleaned input paths
MAY_OUTFLOW_PATH = '//MAY/FLOW-OUT_PEAKCLEAN_2009_2026//1HOUR/CWMS/'
MAY_STOR_PATH    = '//MAY/STOR_PEAKCLEAN_2009_2026//1HOUR/CWMS/'

pathnames = [
    '//MOS/FLOW-IN//1HOUR/CWMS/',
    '//MOS/FLOW-OUT//1HOUR/CWMS/',
    '//MOS/STOR//1HOUR/CWMS/',
    '//MAY/FLOW-IN//1HOUR/CWMS/',
    MAY_OUTFLOW_PATH,
    MAY_STOR_PATH,
] + list(ref_paths.values()) + list(local_flow_paths.values())

startDate = "01Oct2008 00:00:00"
endDate   = "01May2026 00:00:00"

# --- Read ---
DataDict = read_from_dss(dss_file="../../CAS_Unreg_FF/data/obsData.dss", pathnames=pathnames,
                         startDate=startDate, endDate=endDate)

# --- Compute local unregulated flow: Castle Rock - Mayfield - Toutle ---
castle_rock_df = DataDict.get(local_flow_paths['Cowlitz at Castle Rock'])
mayfield_df    = DataDict.get(local_flow_paths['Cowlitz below Mayfield Dam'])
toutle_df      = DataDict.get(ref_paths['Toutle at Tower Rd'])

if all(df is not None and not df.empty for df in [castle_rock_df, mayfield_df, toutle_df]):
    castle   = castle_rock_df.iloc[:, 0].replace(-902, np.nan)
    mayfield = mayfield_df.iloc[:, 0].replace(-902, np.nan)
    toutle   = toutle_df.iloc[:, 0].replace(-902, np.nan)

    local_flow = (castle - mayfield - toutle).clip(lower=0)
    local_path = '/LOCAL/CASTLE_ROCK_LOCAL/FLOW//1HOUR/DERIVED/'
    DataDict[local_path] = local_flow.rename('value').to_frame()
    ref_paths['Castle Rock Local'] = local_path
    print(f"[INFO] Computed Castle Rock local flow: "
          f"{local_flow.notna().sum()} valid values, "
          f"mean={local_flow.mean():.0f} cfs, max={local_flow.max():.0f} cfs")
else:
    print("[WARNING] Could not compute Castle Rock local flow — missing input data.")

# --- Calculate inflow variants ---
# MOS: standard paths
DataDict = recalculate_inflow_variants(DataDict, 'MOS')

# MAY: use cleaned peak paths for storage and outflow
DataDict = recalculate_inflow_variants(DataDict, 'MAY',
                                       stor_path_override=MAY_STOR_PATH,
                                       outflow_path_override=MAY_OUTFLOW_PATH)

# --- Diagnostic ---
for reservoir in ['MOS', 'MAY']:
    orig = DataDict.get(f'//{reservoir}/FLOW-IN//1HOUR/CWMS/')
    if orig is not None:
        print(f"\n{reservoir} FLOW-IN: {orig.index.min()} to {orig.index.max()}, "
              f"{orig['value'].notna().sum()} non-NaN")

# --- Compute combined peaks once, reuse everywhere ---
ref_series_for_combined = {}
for label, path in ref_paths.items():
    df = DataDict.get(path)
    if df is not None and not df.empty:
        ref_series_for_combined[label] = df.iloc[:, 0].replace(-902, np.nan).dropna()

combined_peaks = compute_combined_peaks(ref_series_for_combined, window_hours=36)

# --- Full record inflow comparison plot ---
plot_inflow_interactive(DataDict, ['MOS', 'MAY'])

# --- Peak window plots ---
plot_peak_windows(DataDict, reservoir='MOS', ref_paths=ref_paths,
                  combined_peaks=combined_peaks, window_days=15)
plot_peak_windows(DataDict, reservoir='MAY', ref_paths=ref_paths,
                  combined_peaks=combined_peaks, window_days=15)

# --- Write peak window data to DSS ---
for reservoir in ['MOS', 'MAY']:
    write_peak_windows_to_dss(
        dss_file='../../CAS_Unreg_FF/data/obsData.dss',
        DataDict=DataDict,
        reservoir=reservoir,
        ref_paths=ref_paths,
        combined_peaks=combined_peaks,
        window_days=15
    )

# --- Write full record raw calculated inflow to DSS ---
print("\n[INFO] Writing full record raw inflow to DSS...")
for reservoir in ['MOS', 'MAY']:
    write_full_record_to_dss(
        dss_file='../../CAS_Unreg_FF/data/obsData.dss',
        DataDict=DataDict,
        reservoir=reservoir,
        pathname_suffix='FLOW-IN-CALC-RAW'
    )