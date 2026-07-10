"""
cwms_io.py

Shared helpers for downloading CWMS data via the cwms library (USACE
NWDP data API). Flat functions, no classes. Import what you need:

    from Modules.cwms_io import cwms_download, process_cwms_data

Conventions assumed across projects:
- sites_dict maps {tsid: friendly_name}.
- office defaults to "nwp" (Portland District); pass a different office
  code (e.g. "nws") for other districts.
- Sentinel values -901, -902, and anything below -9000 are treated as gaps.
"""

import pandas as pd
import cwms


def init_cwms_session(office="nwp"):
    """
    Initialize a CWMS API session for a given district office.

    Parameters
    ----------
    office : str
        District office code, e.g. "nwp" (Portland), "nws" (Seattle).

    Returns
    -------
    cwms API object.
    """
    api_root = f"https://wm.{office}.ds.usace.army.mil:8243/nwdp-data/"
    api = cwms.api.init_session(api_root=api_root)
    return api


def cwms_download(sites_dict, start_date, end_date, office="nwp"):
    """
    Download CWMS time series data for a set of tsids.

    Parameters
    ----------
    sites_dict : dict
        {tsid: friendly_name}.
    start_date, end_date : str or datetime-like
        Period to download.
    office : str
        District office code passed to init_cwms_session and get_timeseries.

    Returns
    -------
    dict of {friendly_name: raw dataframe}. Failed/empty sites are omitted.
    """
    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)

    init_cwms_session(office=office)

    results = {}
    for site, name in sites_dict.items():
        try:
            data = cwms.get_timeseries(
                site, office_id=office.upper(), begin=start_date, end=end_date
            ).df
            if data.empty:
                print(f"Downloaded data for {site} is empty.")
            else:
                results[name] = data
        except Exception as e:
            print(f"Failed to download data for {site}: {e}")
    return results


def cwms_download_raw(sites_dict, start_date, end_date, office="nwp", value_col="value"):
    """
    Download CWMS data and return it unresampled/irregular, for records
    that shouldn't be forced onto a regular grid (e.g. event-based data,
    irregular sensor reporting).

    Parameters
    ----------
    sites_dict : dict
        {tsid: friendly_name}.
    start_date, end_date : str or datetime-like
        Period to download.
    office : str
        District office code.
    value_col : str
        Name of the column holding values in the raw CWMS dataframe.

    Returns
    -------
    dict of {friendly_name: pandas.Series}, native timestamps, sentinel
    values masked out and dropped (not filled).
    """
    raw = cwms_download(sites_dict, start_date, end_date, office=office)
    results = {}
    for name, df in raw.items():
        series = pd.to_numeric(df[value_col], errors="coerce")
        series = series.mask(series.isin([-901, -902]))
        series = series.mask(series < -9000)
        results[name] = series.dropna()
    return results


def process_cwms_data(df, value_col, target_interval="1HOUR"):
    """
    Clean and resample a raw CWMS dataframe to a regular interval series.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw dataframe returned by cwms_download.
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

    # mask sentinel / bad values (CWMS downloads can include odd formats)
    series = series.mask(series.isin([-901, -902]))
    series = series.mask(series < -9000)

    freq = "H" if "HOUR" in target_interval.upper() else "D"
    resampled = series.resample(freq).mean()

    resampled = resampled.fillna(-902)

    return resampled
