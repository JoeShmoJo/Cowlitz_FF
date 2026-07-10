"""
dss_io.py

Shared read/write helpers for HEC-DSS files via pydsstools.
Flat functions, no classes. Import what you need:

    from Modules.dss_io import read_ts, write_regular_ts, read_dss_paths

Conventions assumed across projects:
- Pathnames follow A/B/C/D/E/F CWMS-style parts.
- Units are inferred from pathname content (FLOW, ELEV, STOR) unless
  explicitly passed in.
- Missing/sentinel values (-901, -902, below -9000) are treated as gaps.
"""

from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer
import pandas as pd
import numpy as np


def infer_units(pathname):
    """Infer DSS units from pathname content (Part C or full string)."""
    p = pathname.upper()
    if "FLOW" in p:
        return "cfs"
    elif "ELEV" in p:
        return "feet"
    elif "STOR" in p:
        return "ac-ft"
    else:
        return "UNDEFINED"


def infer_type(pathname):
    """Infer DSS data type (PER-AVE vs INST-VAL) from pathname content."""
    p = pathname.upper()
    if "FLOW" in p:
        return "PER-AVE"
    else:
        return "INST-VAL"


def read_ts(dss_path, pathname):
    """
    Read a time series from a DSS file into a pandas Series. Works for
    both regular and irregular records -- pydsstools' read_ts returns
    whatever the record actually is based on the DSS file's own
    metadata, so no separate call is needed for irregular data.

    Parameters
    ----------
    dss_path : str
        Path to the .dss file.
    pathname : str
        Full CWMS-style pathname (A/B/C/D/E/F).

    Returns
    -------
    pandas.Series indexed by datetime, or None if the read fails.
    """
    try:
        with HecDss.Open(dss_path) as fid:
            ts = fid.read_ts(pathname)
            dates = pd.to_datetime(ts.pytimes)
            values = np.array(ts.values, dtype=float)
            series = pd.Series(values, index=dates, name=pathname)
        return series
    except Exception as e:
        print(f"Failed to read {pathname} from {dss_path}: {e}")
        return None


# kept as an alias -- read_ts handles both regular and irregular records
read_regular_ts = read_ts


def read_dss_paths(dss_path, pathnames):
    """
    Read multiple pathnames from one DSS file.

    Parameters
    ----------
    dss_path : str
        Path to the .dss file.
    pathnames : list of str
        Full CWMS-style pathnames.

    Returns
    -------
    dict of {pathname: pandas.Series}. Failed reads are omitted.
    """
    results = {}
    with HecDss.Open(dss_path) as fid:
        for pathname in pathnames:
            try:
                ts = fid.read_ts(pathname)
                dates = pd.to_datetime(ts.pytimes)
                values = np.array(ts.values, dtype=float)
                results[pathname] = pd.Series(values, index=dates, name=pathname)
            except Exception as e:
                print(f"Failed to read {pathname}: {e}")
    return results


def write_regular_ts(dss_path, pathname, series, interval, units=None, dtype=None):
    """
    Write a regular-interval pandas Series to a DSS file.

    Parameters
    ----------
    dss_path : str
        Path to the .dss file (created if it doesn't exist).
    pathname : str
        Full CWMS-style pathname (A/B/C/D/E/F).
    series : pandas.Series
        Datetime-indexed series of values, regular interval.
    interval : str
        DSS interval string, e.g. "1Hour", "1Day".
    units : str, optional
        DSS units string. Inferred from pathname if not given.
    dtype : str, optional
        DSS data type ("PER-AVE" or "INST-VAL"). Inferred from pathname
        if not given.
    """
    if units is None:
        units = infer_units(pathname)
    if dtype is None:
        dtype = infer_type(pathname)

    tsc = TimeSeriesContainer()
    tsc.pathname = pathname
    tsc.startDateTime = series.index[0].strftime("%d%b%Y %H:%M")
    tsc.numberValues = len(series)
    tsc.units = units
    tsc.type = dtype
    tsc.interval = interval
    tsc.values = series.values.astype(float)

    with HecDss.Open(dss_path) as fid:
        fid.put_ts(tsc)


def write_irregular_ts(dss_path, pathname, series, units=None, dtype="INST-VAL"):
    """
    Write an irregular-interval pandas Series (e.g. peaks, EMA input)
    to a DSS file.

    Parameters
    ----------
    dss_path : str
        Path to the .dss file.
    pathname : str
        Full CWMS-style pathname (A/B/C/D/E/F).
    series : pandas.Series
        Datetime-indexed series of values, irregular timestamps.
    units : str, optional
        DSS units string. Inferred from pathname if not given.
    dtype : str, optional
        DSS data type. Defaults to "INST-VAL" for irregular series.
    """
    if units is None:
        units = infer_units(pathname)

    tsc = TimeSeriesContainer()
    tsc.pathname = pathname
    tsc.interval = 0
    tsc.numberValues = len(series)
    tsc.units = units
    tsc.type = dtype
    tsc.times = list(series.index.to_pydatetime())
    tsc.values = series.values.astype(float)

    with HecDss.Open(dss_path) as fid:
        fid.put_ts(tsc)


def parse_pathname(pathname):
    """
    Split a full CWMS-style DSS pathname into its A/B/C/D/E/F parts.

    Parameters
    ----------
    pathname : str
        e.g. "//MAY/FLOW-OUT_PEAKCLEAN_2009_2026/01JAN2009/1HOUR//"

    Returns
    -------
    dict with keys a, b, c, d, e, f
    """
    parts = pathname.strip("/").split("/") if pathname.startswith("/") else pathname.split("/")
    # pathname format is /A/B/C/D/E/F/ -> split("/") on the full string
    # (including leading/trailing slashes) gives ['', A, B, C, D, E, F, '']
    raw = pathname.split("/")
    if len(raw) < 8:
        raise ValueError(f"Pathname does not look like a valid A-F path: {pathname}")
    return {
        "a": raw[1],
        "b": raw[2],
        "c": raw[3],
        "d": raw[4],
        "e": raw[5],
        "f": raw[6],
    }


def build_pathname(a="", b="", c="", d="", e="", f=""):
    """Build a full CWMS-style DSS pathname from individual A-F parts."""
    return f"/{a}/{b}/{c}/{d}/{e}/{f}/"
