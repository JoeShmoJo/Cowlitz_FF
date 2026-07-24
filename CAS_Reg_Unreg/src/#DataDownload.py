# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:00:45 2024
# Change Log
- 04Oct2024
Removed Timestep as a variable and now resample to hourly or daily based on defined ResSim path in the process USGS and process CWMS functions 
Added df = df.apply(pd.to_numeric,errors='coerce') to process CWMS function to catch odd format values downloading from CWMS  

- 05Dec2025
Replaced wincertstore-based SSL handling with a combined CA bundle built
from certifi + Windows ROOT store using ssl.enum_certificates (Py3.9+).
Falls back to certifi only if Windows store access is unavailable.

- May2026
Removed interactive date prompt; dates now defined at top of script.
Added parquet output option.
Fixed all_data build order so parquet gets processed data.

@author: g2encjer
"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import dataretrieval.nwis as nwis
import cwms
from pydsstools.heclib.dss import HecDss

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# NOTE FOR REVIEWERS: this script DOWNLOADS current data and writes DSS/parquet
# into ../../CAS_Unreg_FF/data/obsData.dss (as configured). Re-running regenerates those
# records with freshly downloaded data. NWIS service is chosen per record:
# ResSim paths with a 1DAY E-part use daily values (dv); others use iv.

from pydsstools.core import TimeSeriesContainer

import ssl
import certifi

# =============================================================================
# USER SETTINGS — edit these
# =============================================================================
startDate = "1900-10-01"
endDate   = "2026-10-01"

DSS_OUT_NAME = "../../CAS_Unreg_FF/data/obsData.dss"
PARQUET_OUT_NAME = "../../CAS_Unreg_FF/data/obsData.parquet"


# Output format: 'parquet', 'dss', or 'both'
output_format = 'dss'
# =============================================================================

# Base folder where the script/EXE lives
if getattr(sys, "frozen", False):
    BASE_DIR = Path(sys.argv[0]).resolve().parent
else:
    BASE_DIR = Path(__file__).resolve().parent

RequiredRecordsDictPath = BASE_DIR / '../data/MOS_ELEV.csv'
ObsDataWrite = str(BASE_DIR / DSS_OUT_NAME)
ObsDataParquet = str(BASE_DIR / PARQUET_OUT_NAME)

print(f"\nUsing StartDate = {startDate}, EndDate = {endDate}")
print(f"Output format: {output_format}\n")

# --- SSL Certificate Setup ---
pem_path = os.path.join(tempfile.gettempdir(), "corp_plus_certifi.pem")


                             # debug


def build_windows_ca_bundle(target_pem: str) -> str:
    base_bundle = certifi.where()
    with open(base_bundle, "rb") as src, open(target_pem, "wb") as dst:
        dst.write(src.read())
        try:
            for cert_tuple in ssl.enum_certificates("ROOT"):
                der_bytes = cert_tuple[0]
                pem_str = ssl.DER_cert_to_PEM_cert(der_bytes)
                dst.write(pem_str.encode("ascii"))
        except AttributeError:
            print("[WARNING] ssl.enum_certificates not available; using certifi only.")
        except Exception as e:
            print(f"[WARNING] Error reading Windows ROOT store: {e}")
    return target_pem

if not os.path.exists(pem_path):
    try:
        bundle_path = build_windows_ca_bundle(pem_path)
        print(f"[INFO] Built combined CA bundle: {bundle_path}")
    except Exception as e:
        print(f"[WARNING] Failed to build combined CA bundle: {e}")
        bundle_path = certifi.where()
else:
    bundle_path = pem_path

os.environ["REQUESTS_CA_BUNDLE"] = bundle_path
print(f"[INFO] Using CA bundle: {bundle_path}")

# --- End SSL Setup ---


# --- Functions ---

def NWIS_dl(sites_dict, service, startDate, endDate, parameterCD):
    NWIS = {}
    for site, name in sites_dict.items():
        try:
            data = nwis.get_record(
                sites=site,
                service=service,
                start=startDate,
                end=endDate,
                parameterCd=parameterCD
            )
            if data.empty:
                print(f"Downloaded data for {site} is empty.")
            else:
                NWIS[name] = data
                print(f"Successfully downloaded data for {site}.")
        except Exception as e:
            print(f"Failed to download data for {site}: {e}")
    return NWIS


def CWMS_Download(sites_dict, StartDate, EndDate, office='nws'):
    StartDate = pd.to_datetime(StartDate)
    EndDate = pd.to_datetime(EndDate)
    apiRoot = "https://wm." + office + ".ds.usace.army.mil:8243/nwdp-data/"
    api = cwms.api.init_session(api_root=apiRoot)
    CWMS_data = {}
    for site, name in sites_dict.items():
        try:
            data = cwms.get_timeseries(
                site,
                office_id='NWDP',
                begin=StartDate,
                end=EndDate
            ).df
            if data.empty:
                print(f"Downloaded data for {site} is empty.")
            else:
                CWMS_data[name] = data
                print(f"Successfully downloaded data for {site}.")
        except Exception as e:
            print(f"Failed to download data for {site}: {e}")
    return CWMS_data


def process_usgs_data(DataDict):
    results = []
    for df_name, df in DataDict.items():
        if not df.empty and df.shape[1] > 0:
            # Select the flow/elev column by parameter code, not position
            data_cols = [c for c in df.columns if str(c).startswith('00060') or str(c).startswith('62614') or str(c).startswith('00062')]
            data_cols = [c for c in data_cols if not c.endswith('_cd')]  # drop qualifier columns
            if not data_cols:
                print(f"[WARNING] {df_name} has no recognizable data column; skipping.")
                continue
            df = df[[data_cols[0]]].copy()
            df = df.apply(pd.to_numeric, errors='coerce')
            df = df.dropna()
            if df.empty:
                print(f"[WARNING] {df_name} has no numeric data; skipping.")
                continue
            df[df < -9000] = np.nan
            df[df == -902] = np.nan
            df[df == -901] = np.nan
            df = df.dropna()
            first_timestamp = df.index.min().strftime('%Y-%m-%d %H:%M')
            last_timestamp = df.index.max().strftime('%Y-%m-%d %H:%M')
            time_diffs = df.index.to_series().diff().dropna()
            max_gap = time_diffs.max()
            max_gap_hours = max_gap.total_seconds() / 3600.0
            results.append({
                'DataFrame': df_name,
                'First Timestamp': first_timestamp,
                'Last Timestamp': last_timestamp,
                'Max Gap': max_gap,
                'Max Gap Hours': max_gap_hours
            })
            if '1HOUR' in df_name:
                t = 'h'
            elif '1DAY' in df_name:
                t = 'D'
            else:
                print('timestep of ResSim path not 1HOUR or 1DAY')
                continue
            df = df.resample(t).mean()
            df = df.fillna(-902)
            DataDict[df_name] = df
        else:
            print(f"DataFrame {df_name} is either empty or does not have any columns.")
    return pd.DataFrame(results)

# def process_usgs_data(DataDict):
#     results = []
#     for df_name, df in DataDict.items():
#         print(f"\n--- {df_name} ---")
#         print(f"columns: {df.columns.tolist()}")
#         print(f"dtypes:\n{df.dtypes}")
#         print(df.head(2))
#     return pd.DataFrame(results)  

def process_cwms_data(DataDict):
    results = []
    for df_name, df in DataDict.items():
        if 'date-time' in df.columns:
            df = df.set_index('date-time')
            DataDict[df_name] = df
        if not df.empty and df.shape[1] > 0:
            df = df.iloc[:, [0]].copy()
            df = df.apply(pd.to_numeric, errors='coerce')
            df = df.dropna()
            if df.empty:
                print(f"[WARNING] {df_name} has no numeric data; skipping.")
                continue
            df[df < -9000] = np.nan
            df[df == -902] = np.nan
            df[df == -901] = np.nan
            df = df.dropna()
            first_timestamp = df.index.min().strftime('%Y-%m-%d %H:%M')
            last_timestamp = df.index.max().strftime('%Y-%m-%d %H:%M')
            time_diffs = df.index.to_series().diff().dropna()
            max_gap = time_diffs.max()
            max_gap_hours = max_gap.total_seconds() / 3600.0
            results.append({
                'DataFrame': df_name,
                'First Timestamp': first_timestamp,
                'Last Timestamp': last_timestamp,
                'Max Gap': max_gap,
                'Max Gap Hours': max_gap_hours
            })
            if '1HOUR' in df_name:
                t = 'h'
            elif '1DAY' in df_name:
                t = 'D'
            else:
                print('timestep of ResSim path not 1HOUR or 1DAY')
                continue
            df = df.resample(t).mean()
            df = df.fillna(-902)
            DataDict[df_name] = df
        else:
            print(f"DataFrame {df_name} is either empty or does not have any columns.")
    return pd.DataFrame(results)


def write_to_dss(dss_file, DataDict):
    for pathname, df in DataDict.items():
        if df.empty:
            print(f"[WARNING] {pathname} is empty; skipping DSS write.")
            continue
        tsc = TimeSeriesContainer()
        tsc.pathname = pathname
        tsc.startDateTime = str(df.index[0])
        tsc.numberValues = df.shape[0]
        tsc.values = df.iloc[:, 0].copy().to_numpy()
        tsc.interval = 1
        if "ELEV" in pathname.upper():
            tsc.units = "FEET"
        elif "FLOW" in pathname.upper():
            tsc.units = "CFS"
        elif "STOR" in pathname.upper():
            tsc.units = "ACRE-FT"
        else:
            tsc.units = 'Unknown'
            print(f'[WARNING] Path {pathname} is not Flow, Elev, or Stor!')
        tsc.type = "INST-VAL"
        with HecDss.Open(dss_file, version=6) as fid:
            fid.put_ts(tsc)


def write_to_parquet(parquet_file, DataDict):
    frames = []
    for pathname, df in DataDict.items():
        if df.empty:
            print(f"[WARNING] {pathname} is empty; skipping parquet write.")
            continue
        df_out = df.iloc[:, [0]].copy()
        df_out.columns = ['value']
        df_out['pathname'] = pathname
        frames.append(df_out)

    if frames:
        combined = pd.concat(frames)
        combined.to_parquet(parquet_file)
        print(f"[INFO] Wrote parquet: {parquet_file}  ({len(frames)} series)")
    else:
        print("[WARNING] No data to write to parquet.")


def read_from_parquet(parquet_file):
    combined = pd.read_parquet(parquet_file)
    DataDict = {}
    for pathname, group in combined.groupby('pathname'):
        DataDict[pathname] = group[['value']].copy()
    return DataDict


# --- Load required records ---
RequiredRecordsDict = pd.read_csv(RequiredRecordsDictPath, dtype=str)

USGS_df = RequiredRecordsDict[RequiredRecordsDict['Source'] == 'USGS']

USGS_Elev_df = USGS_df[USGS_df['ResSimPath'].str.contains('ELEV', na=False)]
USGS_Elev_dict = dict(zip(USGS_Elev_df['Download_Key'], USGS_Elev_df['ResSimPath']))

USGS_Flow_df = USGS_df[USGS_df['ResSimPath'].str.contains('FLOW', na=False)]
USGS_Flow_dict = dict(zip(USGS_Flow_df['Download_Key'], USGS_Flow_df['ResSimPath']))

CWMS_df = RequiredRecordsDict[RequiredRecordsDict['Source'] == 'CWMS']
CWMS_dict = dict(zip(CWMS_df['Download_Key'], CWMS_df['ResSimPath']))

# --- Download ---
# Daily-value (dv) service for records whose ResSim path is 1DAY; iv otherwise.
USGS_Elev_dv = {k: v for k, v in USGS_Elev_dict.items() if '1DAY' in v}
USGS_Elev_iv = {k: v for k, v in USGS_Elev_dict.items() if '1DAY' not in v}
USGS_Flow_dv = {k: v for k, v in USGS_Flow_dict.items() if '1DAY' in v}
USGS_Flow_iv = {k: v for k, v in USGS_Flow_dict.items() if '1DAY' not in v}

USGS_Elev_Data_Dict = {}
if USGS_Elev_iv:
    USGS_Elev_Data_Dict.update(NWIS_dl(
        sites_dict=USGS_Elev_iv, service='iv',
        startDate=startDate, endDate=endDate, parameterCD='62614'))
if USGS_Elev_dv:
    got = NWIS_dl(sites_dict=USGS_Elev_dv, service='dv',
                  startDate=startDate, endDate=endDate, parameterCD='62614')
    # Legacy fallback: some historical reservoir records are stored under
    # 00062 (elevation above gage datum) instead of 62614 (NAVD88).
    missing = {k: v for k, v in USGS_Elev_dv.items() if v not in got}
    if missing:
        print(f"[INFO] retrying {len(missing)} elevation dv site(s) with "
              f"parameter 00062 (gage-datum elevations -- verify datum "
              f"before use)")
        got.update(NWIS_dl(sites_dict=missing, service='dv',
                           startDate=startDate, endDate=endDate,
                           parameterCD='00062'))
    USGS_Elev_Data_Dict.update(got)

USGS_Flow_Data_Dict = {}
if USGS_Flow_iv:
    USGS_Flow_Data_Dict.update(NWIS_dl(
        sites_dict=USGS_Flow_iv, service='iv',
        startDate=startDate, endDate=endDate, parameterCD='00060'))
if USGS_Flow_dv:
    USGS_Flow_Data_Dict.update(NWIS_dl(
        sites_dict=USGS_Flow_dv, service='dv',
        startDate=startDate, endDate=endDate, parameterCD='00060'))
CWMS_Data_Dict = CWMS_Download(
    sites_dict=CWMS_dict,
    StartDate=startDate,
    EndDate=endDate
)

# --- Process ---
CWMS_Summary_Stats = process_cwms_data(CWMS_Data_Dict)
USGS_Flow_Summary_Stats = process_usgs_data(USGS_Flow_Data_Dict)
USGS_Elev_Summary_Stats = process_usgs_data(USGS_Elev_Data_Dict)

# Build all_data AFTER processing so it contains cleaned, processed data
all_data = {**USGS_Flow_Data_Dict, **USGS_Elev_Data_Dict, **CWMS_Data_Dict}

# --- Summary stats ---
# --- Summary stats ---
summary_frames = [df for df in [CWMS_Summary_Stats, USGS_Flow_Summary_Stats, USGS_Elev_Summary_Stats] if not df.empty]

if summary_frames:
    Combined_Summary_Stats = pd.concat(summary_frames, ignore_index=True)
    Combined_Summary_Stats['Max Gap Hours'] = Combined_Summary_Stats['Max Gap Hours'].astype(float).round(2)
    Combined_Summary_Stats.sort_values(by='Max Gap Hours', ascending=False, inplace=True)
    Combined_Summary_Stats.reset_index(drop=True, inplace=True)
    Combined_Summary_Stats.to_csv('../output/Combined_Summary_Stats.csv', index=None)
else:
    print("[WARNING] No summary stats to write — all datasets were empty or skipped.")
# --- Write output ---
if output_format in ('dss', 'both'):
    write_to_dss(dss_file=ObsDataWrite, DataDict=USGS_Flow_Data_Dict)
    write_to_dss(dss_file=ObsDataWrite, DataDict=USGS_Elev_Data_Dict)
    write_to_dss(dss_file=ObsDataWrite, DataDict=CWMS_Data_Dict)

if output_format in ('parquet', 'both'):
    write_to_parquet(parquet_file=ObsDataParquet, DataDict=all_data)

import plotly.graph_objects as go
import webbrowser
import os
def plot_data_coverage(all_data):
    fig = go.Figure()

    # Get global date range across all series
    global_min = min(df.index.min() for df in all_data.values() if not df.empty)
    global_max = max(df.index.max() for df in all_data.values() if not df.empty)

    for i, (pathname, df) in enumerate(all_data.items()):
        parts = pathname.strip('/').split('/')
        label = ' / '.join([p for p in parts if p])

        # Always draw the full red bar
        fig.add_trace(go.Scatter(
            x=[global_min, global_max],
            y=[label, label],
            mode='lines',
            line=dict(width=8, color='red'),
            name=label,
            showlegend=False,
            hoverinfo='skip',
        ))

        # Skip blue overlay if no data
        if df.empty:
            continue

        valid = df.iloc[:, 0].copy()
        valid = pd.to_numeric(valid, errors='coerce')
        valid = valid.replace(-902, np.nan)
        valid_mask = valid.notna()

        if valid_mask.sum() == 0:
            continue

        # Valid data blocks as blue overlay
        changes = valid_mask.astype(int).diff().fillna(0)
        starts = df.index[changes == 1].tolist()
        ends   = df.index[changes == -1].tolist()

        if valid_mask.iloc[0]:
            starts = [df.index[0]] + starts
        if valid_mask.iloc[-1]:
            ends = ends + [df.index[-1]]

        for start, end in zip(starts, ends):
            fig.add_trace(go.Scatter(
                x=[start, end],
                y=[label, label],
                mode='lines',
                line=dict(width=8, color='steelblue'),
                name=label,
                showlegend=False,
                hovertemplate=f'{label}<br>%{{x}}<extra></extra>',
            ))

    # Manual legend
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='lines',
                             line=dict(width=8, color='steelblue'),
                             name='Data present', showlegend=True))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='lines',
                             line=dict(width=8, color='red'),
                             name='Missing data', showlegend=True))

    fig.update_layout(
        title='Data Coverage by Site',
        xaxis_title='Date',
        yaxis_title='Site',
        height=max(400, len(all_data) * 35),
        hovermode='closest',
        margin=dict(l=400),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
    )

    output_path = os.path.abspath('../diagnostics/data_coverage.html')
    fig.write_html(output_path)
    print(f"Saved: {output_path}")
    webbrowser.open(f'file:///{output_path}')


plot_data_coverage(all_data)