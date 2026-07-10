"""
nwis_io.py

Shared helpers for downloading USGS NWIS data via dataretrieval.
Flat functions, no classes. Import what you need:

    from Modules.nwis_io import nwis_download, process_usgs_data

Conventions assumed across projects:
- sites_dict maps {site_number: friendly_name}.
- parameterCd 00060 = discharge (cfs), 62614 = lake elevation (feet).
- Sentinel values -901, -902, and anything below -9000 are treated as gaps.
"""

import pandas as pd
import dataretrieval.nwis as nwis


def nwis_download(sites_dict, service, start_date, end_date, parameter_cd):
    """
    Download NWIS data for a set of sites.

    Parameters
    ----------
    sites_dict : dict
        {site_number: friendly_name}.
    service : str
        NWIS service, e.g. "iv" (instantaneous values) or "dv" (daily values).
    start_date, end_date : str
        Date strings, e.g. "2023-12-25".
    parameter_cd : str
        USGS parameter code, e.g. "00060" for discharge.

    Returns
    -------
    dict of {friendly_name: raw dataframe}. Failed/empty sites are omitted.
    """
    results = {}
    for site, name in sites_dict.items():
        try:
            data = nwis.get_record(
                sites=site, service=service, start=start_date, end=end_date,
                parameterCd=parameter_cd
            )
            if data.empty:
                print(f"Downloaded data for {site} is empty.")
            else:
                results[name] = data
        except Exception as e:
            print(f"Failed to download data for {site}: {e}")
    return results


def nwis_download_peaks(sites_dict):
    """
    Download annual peak flow records for a set of USGS sites. Peaks are
    inherently irregular (one value per year, on the actual peak date),
    so no resampling is applied.

    Parameters
    ----------
    sites_dict : dict
        {site_number: friendly_name}.

    Returns
    -------
    dict of {friendly_name: raw dataframe}, indexed by peak date, with
    NWIS's native columns (peak_va, peak_dt, etc.) intact. Failed/empty
    sites are omitted.
    """
    results = {}
    for site, name in sites_dict.items():
        try:
            data = nwis.get_record(sites=site, service="peaks")
            if data.empty:
                print(f"Downloaded peak data for {site} is empty.")
            else:
                results[name] = data
        except Exception as e:
            print(f"Failed to download peak data for {site}: {e}")
    return results


def process_usgs_peaks(df, value_col="peak_va"):
    """
    Clean a raw NWIS peaks dataframe into an irregular datetime-indexed
    pandas Series, with no resampling (peaks keep their native dates).

    Parameters
    ----------
    df : pandas.DataFrame
        Raw dataframe returned by nwis_download_peaks.
    value_col : str
        Name of the column holding peak values (default "peak_va").

    Returns
    -------
    pandas.Series, irregular datetime index, sentinel/bad values masked
    out and dropped (not filled, since these are discrete events).
    """
    series = df[value_col].copy()
    series = pd.to_numeric(series, errors="coerce")
    series = series.mask(series.isin([-901, -902]))
    series = series.mask(series < -9000)
    series = series.dropna()
    return series


def nwis_download_raw(sites_dict, service, start_date, end_date, parameter_cd):
    """
    Download NWIS data and return it unresampled, for cases where you
    want the native (often irregular / sub-hourly) reporting interval
    rather than a fixed hourly/daily grid.

    Parameters
    ----------
    sites_dict : dict
        {site_number: friendly_name}.
    service : str
        NWIS service, e.g. "iv" (instantaneous values).
    start_date, end_date : str
        Date strings, e.g. "2023-12-25".
    parameter_cd : str
        USGS parameter code, e.g. "00060" for discharge.

    Returns
    -------
    dict of {friendly_name: pandas.Series}, native timestamps, sentinel
    values masked out and dropped (not filled).
    """
    raw = nwis_download(sites_dict, service, start_date, end_date, parameter_cd)
    results = {}
    for name, df in raw.items():
        value_col = df.columns[0]
        series = pd.to_numeric(df[value_col], errors="coerce")
        series = series.mask(series.isin([-901, -902]))
        series = series.mask(series < -9000)
        results[name] = series.dropna()
    return results


def process_usgs_data(df, value_col, target_interval="1HOUR"):
    """
    Clean and resample a raw NWIS dataframe to a regular interval series.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw dataframe returned by nwis_download.
    value_col : str
        Name of the column holding the values of interest.
    target_interval : str
        "1HOUR" or "1DAY" (case-insensitive); controls resample frequency.

    Returns
    -------
    pandas.Series, datetime-indexed, resampled and gap-filled with -902.
    """
    series = df[value_col].copy()
    series = pd.to_numeric(series, errors="coerce")

    # mask sentinel / bad values
    series = series.mask(series.isin([-901, -902]))
    series = series.mask(series < -9000)

    freq = "H" if "HOUR" in target_interval.upper() else "D"
    resampled = series.resample(freq).mean()

    # fill remaining gaps with sentinel
    resampled = resampled.fillna(-902)

    return resampled


def summarize_record(series, name=None):
    """
    Compute period-of-record and max gap stats for a time series.

    Parameters
    ----------
    series : pandas.Series
        Datetime-indexed series, gaps marked with -902.
    name : str, optional
        Label for the record (defaults to series.name).

    Returns
    -------
    dict with keys: name, start, end, max_gap_hours, n_missing
    """
    valid = series[series != -902]
    if valid.empty:
        return {
            "name": name or series.name,
            "start": None,
            "end": None,
            "max_gap_hours": None,
            "n_missing": len(series),
        }

    gaps = series.index.to_series().diff().dt.total_seconds() / 3600.0
    missing_mask = series == -902
    max_gap = gaps[missing_mask].max() if missing_mask.any() else 0.0

    return {
        "name": name or series.name,
        "start": valid.index.min(),
        "end": valid.index.max(),
        "max_gap_hours": max_gap,
        "n_missing": int(missing_mask.sum()),
    }
