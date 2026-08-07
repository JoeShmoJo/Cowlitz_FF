#Inflow_Volume_Correction.py
# -*- coding: utf-8 -*-
"""
Volume-correct a cleaned/smoothed inflow record to match the raw calculated
inflow volumes over each continuous data chunk.

For each contiguous segment in the CLEANED record, scale all values by:
    scale = sum(RAW over segment) / sum(CLEANED over segment)
and write the corrected series back to DSS as FLOW-IN-CALC-CLEANED-VOLCOR.

Note: RAW calculated inflow can have negative values (non-physical but real
artifacts of the dS calculation) and these are preserved in volume calculations.
Only sentinel values (<= -900) are treated as no-data in the raw record.
All negative values in the CLEANED record are treated as no-data.
"""

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
import webbrowser
import os


# =============================================================================
# USER SETTINGS
# =============================================================================
DSS_FILE   = '../../CAS_Unreg_FF/data/obsData.dss'
RESERVOIRS = ['MOS']
startDate  = "01Oct2008 00:00:00"
endDate    = "01May2026 00:00:00"


LOCAL_VOLUME = "//MAY/FLOW-LOCAL//1HOUR/CWMS/"
PERIODS      = "//MOS/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/"
SHAPE        = "/TILTON RIVER ABOVE BEAR CANYON CREEK NEAR CINEBAR, WA/14236200/FLOW//1HOUR/USGS/"
OUTPUT_SHAPED = "//MAY/FLOW-LOCAL-SHAPED//1HOUR/CWMS/"
# =============================================================================


def read_from_dss(dss_file, pathnames, startDate, endDate):
    DataDict = {}
    with HecDss.Open(dss_file) as fid:
        for pathname in pathnames:
            try:
                ts = fid.read_ts(pathname, window=(startDate, endDate))
                times = pd.to_datetime(ts.pytimes)
                values = ts.values
                df = pd.DataFrame({'value': values}, index=times)
                DataDict[pathname] = df
                print(f"Read: {pathname}")
            except Exception as e:
                print(f"Failed to read {pathname}: {e}")
    return DataDict


def write_to_dss(dss_file, pathname, series):
    """Write a single Series to DSS, preserving gaps as -901 sentinels."""
    if series.dropna().empty:
        print(f"[WARNING] {pathname} has no valid data; skipping write.")
        return

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
    print(f"[INFO] Wrote {n_valid} valid + {len(series) - n_valid} sentinel values: "
          f"{pathname} "
          f"({series.index.min().date()} to {series.index.max().date()})")


def get_continuous_segments(series):
    """
    Return a list of (start, end) timestamp pairs for each contiguous
    block of non-NaN values in the series.
    """
    valid = series.notna()
    segments = []
    in_segment = False
    seg_start = None

    for i, (idx, v) in enumerate(valid.items()):
        if v and not in_segment:
            seg_start = idx
            in_segment = True
        elif not v and in_segment:
            segments.append((seg_start, series.index[i - 1]))
            in_segment = False

    if in_segment:
        segments.append((seg_start, series.index[-1]))

    return segments


def volume_correct_cleaned(raw_df, cleaned_df):
    """
    For each contiguous block of valid data in the cleaned record,
    find the matching date range in the raw record and compute:
        scale = sum(raw over range) / sum(cleaned over range)
    then apply the scale to the cleaned segment.
    """
    raw     = raw_df.iloc[:, 0].copy()
    cleaned = cleaned_df.iloc[:, 0].copy()

    raw[raw <= -900] = np.nan
    cleaned[cleaned <= -900] = np.nan
    cleaned[cleaned < 0]     = np.nan

    print(f"\nRaw:     {raw.notna().sum()} valid values  "
          f"min={raw.min():.1f}  max={raw.max():.1f}")
    print(f"Cleaned: {cleaned.notna().sum()} valid values  "
          f"min={cleaned.min():.1f}  max={cleaned.max():.1f}")

    segments = get_continuous_segments(cleaned)
    print(f"\n[INFO] Found {len(segments)} continuous segments in cleaned record:")
    for seg_start, seg_end in segments:
        n = cleaned[seg_start:seg_end].notna().sum()
        print(f"  {seg_start.date()} → {seg_end.date()}  ({n} hours)")

    corrected = cleaned.copy()
    stats = []

    for seg_start, seg_end in segments:
        cleaned_seg = cleaned[seg_start:seg_end]
        raw_seg     = raw[seg_start:seg_end]

        vol_cleaned = cleaned_seg.sum(skipna=True)
        vol_raw     = raw_seg.sum(skipna=True)

        n_cleaned = cleaned_seg.notna().sum()
        n_raw     = raw_seg.notna().sum()

        if n_cleaned == 0:
            print(f"[WARNING] {seg_start.date()} → {seg_end.date()}: "
                  f"no valid cleaned values; skipping.")
            scale = np.nan
        elif vol_cleaned == 0:
            print(f"[WARNING] {seg_start.date()} → {seg_end.date()}: "
                  f"cleaned volume is zero; skipping.")
            scale = np.nan
        else:
            scale = vol_raw / vol_cleaned
            corrected.loc[seg_start:seg_end] = cleaned_seg * scale

        print(f"  {seg_start.date()} → {seg_end.date()}: "
              f"raw={vol_raw:,.0f}  cleaned={vol_cleaned:,.0f}  "
              f"scale={scale:.4f}" if not np.isnan(scale) else
              f"  {seg_start.date()} → {seg_end.date()}: skipped")

        stats.append({
            'Segment Start':        seg_start,
            'Segment End':          seg_end,
            'N Hours Cleaned':      n_cleaned,
            'N Hours Raw':          n_raw,
            'Vol Raw (cfs-hr)':     round(vol_raw, 0),
            'Vol Cleaned (cfs-hr)': round(vol_cleaned, 0),
            'Scale Factor':         round(scale, 6) if not np.isnan(scale) else np.nan,
        })

    return corrected, raw, cleaned, pd.DataFrame(stats)


def shape_volume_to_hydrograph(periods_series, local_vol_series, shape_series):
    """
    For each contiguous segment in periods_series, take the total volume
    from local_vol_series over that same window and redistribute it using
    the shape of shape_series, scaled so the output integrates to that volume.

    Steps per segment:
      1. Sum local_vol over [seg_start, seg_end]  -> target_volume
      2. Slice shape over [seg_start, seg_end], clip negatives to 0
      3. If shape sum > 0: shaped = shape_slice / sum(shape_slice) * target_volume
         If shape sum == 0: shaped = uniform distribution of target_volume
      4. Assemble all segments into a single output series.

    Returns:
        shaped_out  : pd.Series on the union of all segment timestamps
        stats_df    : pd.DataFrame with per-segment diagnostics
    """
    # Clean sentinel values
    periods  = periods_series.copy()
    periods[periods <= -900] = np.nan

    local_vol = local_vol_series.copy()
    local_vol[local_vol <= -900] = np.nan

    shape = shape_series.copy()
    shape[shape <= -900] = np.nan
    shape[shape < 0]     = 0.0   # shape must be non-negative

    segments = get_continuous_segments(periods)
    print(f"\n[INFO] shape_volume_to_hydrograph: "
          f"{len(segments)} segments in PERIODS record")

    output_pieces = []
    stats = []

    for seg_start, seg_end in segments:
        # Target volume from LOCAL_VOLUME over this window
        local_slice  = local_vol[seg_start:seg_end]
        target_vol   = local_slice.sum(skipna=True)
        n_local_valid = local_slice.notna().sum()

        # Shape slice over same window
        shape_slice  = shape[seg_start:seg_end]
        shape_sum    = shape_slice.sum(skipna=True)
        n_shape_valid = shape_slice.notna().sum()

        if n_shape_valid == 0:
            print(f"[WARNING] {seg_start.date()} → {seg_end.date()}: "
                  f"no valid shape values; filling with uniform distribution.")
            n_hours = len(shape_slice)
            uniform_val = target_vol / n_hours if n_hours > 0 else 0.0
            shaped = pd.Series(uniform_val, index=shape_slice.index)
        elif shape_sum == 0:
            print(f"[WARNING] {seg_start.date()} → {seg_end.date()}: "
                  f"shape sums to zero; filling with uniform distribution.")
            n_hours = n_shape_valid
            uniform_val = target_vol / n_hours if n_hours > 0 else 0.0
            shaped = shape_slice.copy()
            shaped[shaped.notna()] = uniform_val
        else:
            # Normalise shape then scale to target volume
            shaped = (shape_slice / shape_sum) * target_vol

        output_pieces.append(shaped)

        print(f"  {seg_start.date()} → {seg_end.date()}: "
              f"target_vol={target_vol:,.0f} cfs-hr  "
              f"shape_sum={shape_sum:,.0f}  "
              f"shaped_sum={shaped.sum(skipna=True):,.0f}  "
              f"n_local={n_local_valid}  n_shape={n_shape_valid}")

        stats.append({
            'Segment Start':          seg_start,
            'Segment End':            seg_end,
            'N Hours Periods':        len(periods[seg_start:seg_end]),
            'N Hours Local Valid':    n_local_valid,
            'N Hours Shape Valid':    n_shape_valid,
            'Target Vol (cfs-hr)':    round(target_vol, 0),
            'Shape Sum (cfs-hr)':     round(shape_sum, 0),
            'Output Vol (cfs-hr)':    round(shaped.sum(skipna=True), 0),
        })

    if not output_pieces:
        print("[WARNING] No segments produced any output.")
        return pd.Series(dtype=float), pd.DataFrame(stats)

    shaped_out = pd.concat(output_pieces).sort_index()
    # Drop duplicate indices if segments happened to overlap
    shaped_out = shaped_out[~shaped_out.index.duplicated(keep='first')]

    return shaped_out, pd.DataFrame(stats)


def plot_comparison(raw_series, cleaned_series, corrected_series, reservoir):
    """Plot raw, cleaned, and volume-corrected series for visual QC."""
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=raw_series.index, y=raw_series.values,
        name='Raw Calculated',
        line=dict(color='steelblue', width=1),
        opacity=0.8,
    ))
    fig.add_trace(go.Scatter(
        x=cleaned_series.index, y=cleaned_series.values,
        name='Cleaned (input)',
        line=dict(color='darkorange', width=1, dash='dot'),
        opacity=0.9,
    ))
    fig.add_trace(go.Scatter(
        x=corrected_series.index, y=corrected_series.values,
        name='Cleaned Vol-Corrected',
        line=dict(color='green', width=1.5),
        opacity=0.9,
    ))

    fig.update_layout(
        title=f'{reservoir} — Volume Correction: Raw vs Cleaned vs Corrected',
        xaxis_title='Date',
        yaxis_title='Flow (cfs)',
        height=500,
        hovermode='x unified',
        legend=dict(orientation='h', yanchor='bottom', y=1.02,
                    xanchor='right', x=1),
    )

    output_path = os.path.abspath(f'../output/diagnostics/{reservoir}_volcor_comparison.html')
    fig.write_html(output_path)
    print(f"Saved: {output_path}")
    webbrowser.open(f'file:///{output_path}')


def plot_shaped(local_vol_series, shape_series, shaped_series):
    """Plot LOCAL_VOLUME, SHAPE, and the shaped output for visual QC."""
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=local_vol_series.index, y=local_vol_series.values,
        name='Local Volume (input)',
        line=dict(color='steelblue', width=1),
        opacity=0.8,
    ))
    fig.add_trace(go.Scatter(
        x=shape_series.index, y=shape_series.values,
        name='Shape Record (Tilton)',
        line=dict(color='darkorange', width=1, dash='dot'),
        opacity=0.8,
    ))
    fig.add_trace(go.Scatter(
        x=shaped_series.index, y=shaped_series.values,
        name='Shaped Output',
        line=dict(color='green', width=1.5),
        opacity=0.9,
    ))

    fig.update_layout(
        title='MAY Local Flow — Volume Shaped to Tilton Hydrograph',
        xaxis_title='Date',
        yaxis_title='Flow (cfs)',
        height=500,
        hovermode='x unified',
        legend=dict(orientation='h', yanchor='bottom', y=1.02,
                    xanchor='right', x=1),
    )

    output_path = os.path.abspath('../output/diagnostics/MAY_local_shaped.html')
    fig.write_html(output_path)
    print(f"Saved: {output_path}")
    webbrowser.open(f'file:///{output_path}')


# =============================================================================
# MAIN — volume correction loop
# =============================================================================

for reservoir in RESERVOIRS:
    print(f"\n{'='*60}")
    print(f"Processing {reservoir}")
    print(f"{'='*60}")

    raw_path     = f'//{reservoir}/FLOW-IN-CALC-RAW//1HOUR/CWMS/'
    cleaned_path = f'//{reservoir}/FLOW-IN-CALC-CLEANED//1HOUR/CWMS/'
    output_path  = f'//{reservoir}/FLOW-IN-CALC-CLEANED-VOLCOR//1HOUR/CWMS/'

    DataDict = read_from_dss(
        dss_file=DSS_FILE,
        pathnames=[raw_path, cleaned_path],
        startDate=startDate,
        endDate=endDate,
    )

    raw_df     = DataDict.get(raw_path)
    cleaned_df = DataDict.get(cleaned_path)

    if raw_df is None or raw_df.empty:
        print(f"[ERROR] Raw series not found for {reservoir}; skipping.")
        continue
    if cleaned_df is None or cleaned_df.empty:
        print(f"[ERROR] Cleaned series not found for {reservoir}; skipping.")
        continue

    print("\n--- RAW first 5 timestamps ---")
    print(raw_df.index[:5].tolist())
    print("\n--- CLEANED first 5 timestamps ---")
    print(cleaned_df.index[:5].tolist())
    print("\n--- CLEANED first 5 values ---")
    print(cleaned_df.iloc[:5, 0].tolist())
    c = cleaned_df.iloc[:, 0]
    valid = c[(c > -900) & (c > 0)]
    if not valid.empty:
        print(f"\nFirst valid cleaned timestamp: {valid.index[0]}")
        print(f"Last valid cleaned timestamp:  {valid.index[-1]}")
        print(f"Valid cleaned count: {len(valid)}")

    corrected, raw, cleaned, stats_df = volume_correct_cleaned(raw_df, cleaned_df)

    print(f"\nSegment summary for {reservoir}:")
    print(stats_df.to_string(index=False))
    stats_csv = f'../output/diagnostics/{reservoir}_volcor_stats.csv'
    stats_df.to_csv(stats_csv, index=False)
    print(f"Saved: {stats_csv}")

    write_to_dss(DSS_FILE, output_path, corrected)
    plot_comparison(raw, cleaned, corrected, reservoir)


# =============================================================================
# MAIN — hydrograph shaping
# =============================================================================

print(f"\n{'='*60}")
print("Hydrograph shaping: MAY local flow")
print(f"{'='*60}")

shape_DataDict = read_from_dss(
    dss_file=DSS_FILE,
    pathnames=[PERIODS, LOCAL_VOLUME, SHAPE],
    startDate=startDate,
    endDate=endDate,
)

periods_df   = shape_DataDict.get(PERIODS)
local_vol_df = shape_DataDict.get(LOCAL_VOLUME)
shape_df     = shape_DataDict.get(SHAPE)

if periods_df is None or periods_df.empty:
    print(f"[ERROR] PERIODS record not found ({PERIODS}); skipping shaping.")
elif local_vol_df is None or local_vol_df.empty:
    print(f"[ERROR] LOCAL_VOLUME record not found ({LOCAL_VOLUME}); skipping shaping.")
elif shape_df is None or shape_df.empty:
    print(f"[ERROR] SHAPE record not found ({SHAPE}); skipping shaping.")
else:
    periods_series   = periods_df.iloc[:, 0]
    local_vol_series = local_vol_df.iloc[:, 0]
    shape_series     = shape_df.iloc[:, 0]

    shaped_out, shape_stats_df = shape_volume_to_hydrograph(
        periods_series, local_vol_series, shape_series
    )

    print(f"\nShaping segment summary:")
    print(shape_stats_df.to_string(index=False))
    shape_stats_df.to_csv('../output/diagnostics/MAY_local_shaped_stats.csv', index=False)
    print("Saved: MAY_local_shaped_stats.csv")

    write_to_dss(DSS_FILE, OUTPUT_SHAPED, shaped_out)
    plot_shaped(local_vol_series, shape_series, shaped_out)