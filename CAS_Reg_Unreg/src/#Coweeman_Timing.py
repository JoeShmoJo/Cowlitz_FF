#Coweeman_Timing.py
# -*- coding: utf-8 -*-
"""
Coweeman River against Castle Rock: do the peaks arrive together?

WHY THIS MATTERS
    The Coweeman joins the Cowlitz at Kelso, BELOW Castle Rock, so it is part
    of the ungaged-to-Columbia local area rather than part of the basin the
    regulated curve was developed for. Extending the study downstream by
    coincident frequency needs to know whether a Coweeman flood arrives with
    the Castle Rock flood or after it. If the two peak together the local has
    to be added at full value; if the Coweeman leads or lags by a day the
    combination is much weaker, and assuming coincidence would overstate the
    downstream flow.

    This script measures that lag two independent ways, because neither record
    covers the whole period on its own.

PART 1 -- ANNUAL PEAK RECORDS
    USGS annual instantaneous peaks at both sites, matched by water year. This
    reaches back to the start of the Coweeman record but is coarse: the USGS
    peak file often carries a date with no time, especially before the 1990s.
    Pairs are reported separately depending on whether both sides have a clock
    time, and the no-time pairs are only ever used at day resolution.

PART 2 -- SUB-DAILY RECORDS
    The USGS Coweeman gage (14245000) stops in the 1990s, and the USGS
    instantaneous-value service does not reach back that far anyway, so the
    sub-daily comparison uses a different pair of sources over a later period:
    Washington Department of Ecology 15-minute data at 26C075 for the Coweeman,
    against the USGS instantaneous record at Castle Rock. Events are picked off
    the Castle Rock record, and for each one the lag is measured two ways --
    peak-to-peak, and by cross-correlating the two hydrographs over the event
    window. Cross-correlation is the more robust of the two on a flat or
    double-peaked event, where the argmax can jump hours for no physical
    reason.

NETWORK
    Everything is cached to CACHE_DIR on first download, so a re-run is offline
    and the analysis can be edited without hammering either agency. Delete the
    cache to refresh. The three hosts used are waterservices.usgs.gov,
    nwis.waterdata.usgs.gov and apps.ecology.wa.gov.

THE ECOLOGY FILE FORMAT IS NOT ASSUMED
    Only the daily-value URL was known when this was written:
        .../StationData/Prod/26C075/26C075_2019_DSG_DV.txt
    The 15-minute file name is a guess, so ECOLOGY_SUFFIXES holds several
    candidates and the first that returns data is used and reported. The
    column layout is sniffed rather than hard coded -- delimiter, date column
    and flow column are all detected, and if detection fails the script prints
    the first lines of the file it could not read instead of failing silently.
    If none of the candidates work, paste a known-good 15-minute URL into
    ECOLOGY_URL_OVERRIDE.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import io
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
COWEEMAN_USGS_SITE = "14245000"      # Coweeman River near Kelso, WA
CASTLE_ROCK_SITE = "14243000"        # Cowlitz River at Castle Rock, WA
ECOLOGY_SITE = "26C075"              # Ecology Coweeman gage

CACHE_DIR = r"../data/coweeman"
OUT_DIR = r"../output/diagnostics"
PLOT_STEM = os.path.join(OUT_DIR, "coweeman_timing")

# --- sub-daily period --------------------------------------------------------
# Ecology 15-minute coverage. USGS instantaneous values at Castle Rock start
# around 2007, so the usable overlap is the intersection of the two.
ECOLOGY_FIRST_YEAR = 2006
ECOLOGY_LAST_YEAR = 2020
IV_PARAMETER = "00060"               # discharge, cfs

# --- Ecology endpoint --------------------------------------------------------
ECOLOGY_BASE = ("https://apps.ecology.wa.gov/ContinuousFlowAndWQ/StationData/"
                "Prod/%s/%s_%d_DSG_%s.txt")
# Tried in order; the first that returns parseable data wins and is reported.
# "DV" is the daily file whose URL was known -- it is last so that a 15-minute
# file is preferred, but the script still produces something if only DV exists.
ECOLOGY_SUFFIXES = ["15", "UV", "Inst", "15Min", "QC", "DV"]
# A known-good URL with %d where the year goes, if the candidates all miss.
ECOLOGY_URL_OVERRIDE = None

# --- timezone ----------------------------------------------------------------
# USGS instantaneous values come back timezone-aware and shift with daylight
# saving; Ecology publishes in Pacific Standard Time year round. Mixing the two
# raises on comparison, and worse, an uncorrected hour of DST would land right
# in the middle of the lag being measured. Everything is converted to a fixed
# UTC-8 and then made naive, so both sides are on the same clock all year.
LOCAL_STANDARD_OFFSET_HOURS = -8

# --- event selection ---------------------------------------------------------
N_EVENTS = 25                        # largest independent Castle Rock events
EVENT_MIN_SEPARATION_DAYS = 7        # peaks closer than this are one event
EVENT_WINDOW_HOURS = 96              # +/- window the lag is searched over
RESAMPLE = "1h"                      # common step both series are put on
MIN_EVENT_CFS = 15000.0              # ignore small events at Castle Rock
# A window must have this fraction of its hours present in BOTH series.
MIN_WINDOW_COVERAGE = 0.80
# Cross-correlation search range, hours. Wider than any plausible travel-time
# difference so the optimum is interior rather than pinned at the edge.
XCORR_MAX_LAG_HOURS = 48

# --- annual peak matching ----------------------------------------------------
PEAK_MAX_SEPARATION_DAYS = 5         # beyond this the two are different storms

C_COW = "#b7410e"
C_CAS = "#1a4f8a"

# ----------------------------------------------------------------------------


def ensure_dirs():
    for path in (CACHE_DIR, OUT_DIR):
        if path and not os.path.isdir(path):
            os.makedirs(path)


def cache_path(name):
    return os.path.join(CACHE_DIR, name)


def http_text(url, timeout=120):
    """GET returning text, or None with the reason printed."""
    import requests
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            print("      HTTP %d  %s" % (response.status_code, url))
            return None
        return response.text
    except Exception as exc:
        print("      failed: %s" % exc)
        return None


def cached_text(name, url):
    """Download once, then read from disk."""
    path = cache_path(name)
    if os.path.isfile(path):
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    print("   GET %s" % url)
    text = http_text(url)
    if text is None:
        return None
    with open(path, "w", encoding="utf-8", errors="replace") as handle:
        handle.write(text)
    return text


def parse_rdb(text):
    """USGS RDB to a DataFrame.

    RDB is tab separated with '#' comment lines, a header line, and then a
    FORMAT line ('5s', '15s', ...) that is not data. Dropping only the comments
    and taking the next line as data is the classic way to get one row of
    garbage into the top of the frame.
    """
    if not text:
        return None
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 3:
        return None
    header = lines[0].split("\t")
    body = [ln.split("\t") for ln in lines[2:]]
    body = [r for r in body if len(r) == len(header)]
    if not body:
        return None
    return pd.DataFrame(body, columns=header)


def to_naive_local(index):
    """Any DatetimeIndex to naive local standard time.

    Timezone-aware input is converted to the fixed offset and stripped; naive
    input is assumed to already be on that clock.
    """
    index = pd.DatetimeIndex(index)
    if index.tz is not None:
        offset = pd.Timedelta(hours=LOCAL_STANDARD_OFFSET_HOURS)
        index = index.tz_convert("UTC") + offset
        index = index.tz_localize(None)
    return index


# ----------------------------------------------------------------------------
# USGS
# ----------------------------------------------------------------------------

def fetch_usgs_peaks(site):
    """Annual instantaneous peaks: date, time where published, and value."""
    name = "usgs_peaks_%s.rdb" % site
    urls = [
        "https://nwis.waterdata.usgs.gov/nwis/peak?site_no=%s&agency_cd=USGS"
        "&format=rdb" % site,
        "https://waterdata.usgs.gov/nwis/peak?site_no=%s&agency_cd=USGS"
        "&format=rdb" % site,
    ]
    text = None
    for url in urls:
        text = cached_text(name, url)
        if text and "peak_dt" in text:
            break
        text = None
    table = parse_rdb(text)
    if table is None or "peak_dt" not in table.columns:
        print("   peaks %s: UNAVAILABLE" % site)
        return None

    out = pd.DataFrame()
    out["peak_date"] = pd.to_datetime(table["peak_dt"], errors="coerce")
    time_text = (table["peak_tm"].astype(str).str.strip()
                 if "peak_tm" in table.columns else "")
    out["peak_time_text"] = time_text
    out["peak_cfs"] = pd.to_numeric(table["peak_va"], errors="coerce")
    stamp = out["peak_date"].astype(str) + " " + out["peak_time_text"]
    out["peak_stamp"] = pd.to_datetime(stamp, errors="coerce")
    out["has_time"] = out["peak_stamp"].notna()
    # Where no clock time was published, fall back to the date at midnight so
    # the pair can still be compared at day resolution.
    out["peak_stamp"] = out["peak_stamp"].fillna(out["peak_date"])
    out = out.dropna(subset=["peak_date", "peak_cfs"])
    out["WY"] = out["peak_date"].dt.year + (out["peak_date"].dt.month >= 10)
    print("   peaks %s: %d water years (%d with a clock time), WY%d-%d"
          % (site, len(out), int(out["has_time"].sum()),
             out["WY"].min(), out["WY"].max()))
    return out.reset_index(drop=True)


def fetch_usgs_iv(site, start, end):
    """Instantaneous discharge, cached, via dataretrieval or plain RDB."""
    name = "usgs_iv_%s_%s_%s.csv" % (site, start, end)
    path = cache_path(name)
    if os.path.isfile(path):
        frame = pd.read_csv(path, parse_dates=["datetime"], index_col="datetime")
        series = frame.iloc[:, 0].astype(float)
        series.index = to_naive_local(series.index)
        print("   IV %s: cached (%d values)" % (site, len(series)))
        return series.sort_index()

    series = None
    try:
        import dataretrieval.nwis as nwis
        print("   IV %s: dataretrieval %s to %s" % (site, start, end))
        data = nwis.get_record(sites=site, service="iv", start=start, end=end,
                               parameterCd=IV_PARAMETER)
        if data is not None and len(data):
            column = next((c for c in data.columns if c.startswith(IV_PARAMETER)
                           and not c.endswith("_cd")), None)
            if column:
                series = pd.to_numeric(data[column], errors="coerce").dropna()
    except ImportError:
        print("   IV %s: dataretrieval not installed, using the RDB service"
              % site)
    except Exception as exc:
        print("   IV %s: dataretrieval failed (%s), using the RDB service"
              % (site, exc))

    if series is None:
        url = ("https://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=%s"
               "&startDT=%s&endDT=%s&parameterCd=%s"
               % (site, start, end, IV_PARAMETER))
        table = parse_rdb(cached_text("usgs_iv_%s.rdb" % site, url))
        if table is None:
            print("   IV %s: UNAVAILABLE" % site)
            return None
        value_col = next((c for c in table.columns
                          if c.endswith("_" + IV_PARAMETER)), None)
        if value_col is None:
            value_col = next((c for c in table.columns
                              if IV_PARAMETER in c and not c.endswith("_cd")),
                             None)
        if value_col is None or "datetime" not in table.columns:
            print("   IV %s: no discharge column in %s"
                  % (site, list(table.columns)))
            return None
        stamps = pd.to_datetime(table["datetime"], errors="coerce")
        values = pd.to_numeric(table[value_col], errors="coerce")
        series = pd.Series(values.values, index=stamps).dropna()

    series.index = to_naive_local(series.index)
    series = series[~series.index.duplicated(keep="first")].sort_index()
    series.rename("cfs").to_frame().rename_axis("datetime").to_csv(path)
    print("   IV %s: %d values, %s to %s"
          % (site, len(series), series.index.min(), series.index.max()))
    return series


# ----------------------------------------------------------------------------
# Ecology
# ----------------------------------------------------------------------------

def sniff_ecology(text):
    """Find the date and flow columns in an Ecology station file.

    The layout is not documented here, so it is detected: the delimiter is
    whichever of tab / comma / whitespace parses the most columns, the date
    column is the first that parses as datetimes, and the flow column is the
    first numeric one whose name mentions discharge or flow -- or, failing a
    helpful name, the first numeric column after the date.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None
    start = 0
    for i, line in enumerate(lines[:80]):
        if re.search(r"\d{1,4}[/-]\d{1,2}[/-]\d{1,4}", line):
            start = i
            break
    head = "\n".join(lines[max(start - 1, 0):start + 400])

    best = None
    # sep=r"\s+" rather than delim_whitespace, which pandas 3 removed
    # outright -- it would raise TypeError on any whitespace-delimited file.
    for sep, kwargs in ((",", {"sep": ","}), ("\t", {"sep": "\t"}),
                        ("ws", {"sep": r"\s+"})):
        try:
            trial = pd.read_csv(io.StringIO(head), engine="python", **kwargs)
        except Exception:
            continue
        if trial.shape[1] < 2:
            continue
        if best is None or trial.shape[1] > best[1].shape[1]:
            best = (sep, trial)
    if best is None:
        return None
    sep, frame = best

    date_col = None
    for column in frame.columns:
        parsed = pd.to_datetime(frame[column], errors="coerce")
        if parsed.notna().mean() > 0.8:
            date_col = column
            break
    if date_col is None:
        return None

    # A SEPARATE TIME COLUMN has to be found and joined on, or the whole file
    # collapses to one value per day. Every row then carries the same midnight
    # stamp, the duplicate-index drop keeps one, and 15-minute data quietly
    # becomes daily -- which would not error anywhere, it would just make the
    # sub-daily lag meaningless. Caught on a whitespace-delimited test layout
    # with "Date Time Discharge" columns.
    time_col = None
    for column in frame.columns:
        if column == date_col:
            continue
        text = frame[column].astype(str).str.strip()
        looks_like_time = text.str.match(r"^\d{1,2}:\d{2}(:\d{2})?$").mean()
        if looks_like_time > 0.8:
            time_col = column
            break

    exclude = {date_col, time_col}
    numeric = [c for c in frame.columns if c not in exclude
               and pd.to_numeric(frame[c], errors="coerce").notna().mean() > 0.8]
    if not numeric:
        return None
    named = [c for c in numeric
             if re.search(r"disch|flow|cfs|q\b", str(c), re.I)]
    return {"sep": sep, "date_col": date_col, "time_col": time_col,
            "flow_col": named[0] if named else numeric[0],
            "skip": max(start - 1, 0)}


def parse_ecology(text):
    """One Ecology file to a flow Series, or None."""
    spec = sniff_ecology(text)
    if spec is None:
        return None, None
    lines = [ln for ln in text.splitlines() if ln.strip()]
    body = "\n".join(lines[spec["skip"]:])
    kwargs = ({"sep": r"\s+"} if spec["sep"] == "ws"
              else {"sep": spec["sep"]})
    try:
        frame = pd.read_csv(io.StringIO(body), engine="python", **kwargs)
    except Exception:
        return None, None
    text = frame[spec["date_col"]].astype(str).str.strip()
    if spec.get("time_col"):
        text = text + " " + frame[spec["time_col"]].astype(str).str.strip()
    stamps = pd.to_datetime(text, errors="coerce")
    values = pd.to_numeric(frame[spec["flow_col"]], errors="coerce")
    series = pd.Series(values.values, index=stamps).dropna()
    dropped = int(series.index.duplicated(keep="first").sum())
    series = series[~series.index.duplicated(keep="first")].sort_index()
    if dropped > 0.5 * (len(series) + dropped):
        print("      WARNING: %d of %d stamps were duplicates. If this file is "
              "sub-daily," % (dropped, len(series) + dropped))
        print("      its time column was not recognised and the record has "
              "collapsed to daily.")
    return (series if len(series) else None), spec


def fetch_ecology(site, first_year, last_year):
    """Ecology sub-daily flow across the requested years.

    Each year is tried against ECOLOGY_SUFFIXES until one parses. The suffix
    that worked is reported, because it is a guess -- see the docstring.
    """
    pieces, used, failed = [], {}, []
    for year in range(first_year, last_year + 1):
        got = None
        candidates = ([ECOLOGY_URL_OVERRIDE % year] if ECOLOGY_URL_OVERRIDE
                      else [ECOLOGY_BASE % (site, site, year, s)
                            for s in ECOLOGY_SUFFIXES])
        for url in candidates:
            tag = url.rsplit("_", 1)[-1].replace(".txt", "")
            text = cached_text("ecology_%s_%d_%s.txt" % (site, year, tag), url)
            if not text or len(text) < 200:
                continue
            series, spec = parse_ecology(text)
            if series is not None and len(series) > 24:
                got = (series, tag, spec)
                break
        if got is None:
            failed.append(year)
            continue
        series, tag, spec = got
        pieces.append(series)
        used[year] = (tag, len(series), spec["date_col"], spec["flow_col"])

    if not pieces:
        print("   Ecology %s: NOTHING PARSED for %d-%d"
              % (site, first_year, last_year))
        print("   Tried suffixes: %s" % ", ".join(ECOLOGY_SUFFIXES))
        print("   Set ECOLOGY_URL_OVERRIDE to a known-good URL with %d for the "
              "year.")
        probe = cache_path("ecology_%s_%d_%s.txt"
                           % (site, first_year, ECOLOGY_SUFFIXES[0]))
        if os.path.isfile(probe):
            with open(probe, "r", errors="replace") as handle:
                print("   First lines of what came back:")
                for line in handle.read().splitlines()[:8]:
                    print("      %s" % line[:110])
        return None

    combined = pd.concat(pieces).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]
    combined.index = to_naive_local(combined.index)
    tags = sorted({v[0] for v in used.values()})
    step = pd.Series(combined.index).diff().median()
    print("   Ecology %s: %d values, %s to %s (file type %s, step %s)"
          % (site, len(combined), combined.index.min().date(),
             combined.index.max().date(), "/".join(tags), step))
    if failed:
        print("      no data for: %s" % ", ".join(str(y) for y in failed))
    return combined


# ----------------------------------------------------------------------------
# Part 1 -- annual peaks
# ----------------------------------------------------------------------------

def compare_annual_peaks(cow, cas):
    """Match the two peak records by water year and time the difference."""
    if cow is None or cas is None:
        return None
    merged = cow.merge(cas, on="WY", suffixes=("_cow", "_cas"))
    if not len(merged):
        return None
    merged["lag_days"] = (merged["peak_date_cow"]
                          - merged["peak_date_cas"]).dt.total_seconds() / 86400.0
    both_timed = merged["has_time_cow"] & merged["has_time_cas"]
    merged["both_have_time"] = both_timed
    merged["lag_hours"] = np.where(
        both_timed,
        (merged["peak_stamp_cow"]
         - merged["peak_stamp_cas"]).dt.total_seconds() / 3600.0,
        np.nan)
    merged["same_storm"] = merged["lag_days"].abs() <= PEAK_MAX_SEPARATION_DAYS
    return merged


def report_annual(merged):
    print("\n" + "=" * 78)
    print("PART 1 -- ANNUAL PEAK RECORDS")
    print("=" * 78)
    if merged is None or not len(merged):
        print("   no overlapping water years")
        return
    same = merged[merged["same_storm"]]
    print("   %d water years in common (WY%d-%d); %d within %d days, so the "
          "same storm"
          % (len(merged), merged["WY"].min(), merged["WY"].max(), len(same),
             PEAK_MAX_SEPARATION_DAYS))
    if len(same):
        print("   lag in DAYS (Coweeman minus Castle Rock), same-storm pairs:")
        print("      median %+.2f   mean %+.2f   IQR %+.2f to %+.2f"
              % (same["lag_days"].median(), same["lag_days"].mean(),
                 same["lag_days"].quantile(0.25),
                 same["lag_days"].quantile(0.75)))
        print("      same calendar day: %d of %d (%.0f%%)"
              % (int((same["lag_days"] == 0).sum()), len(same),
                 100.0 * (same["lag_days"] == 0).mean()))
    timed = merged[merged["both_have_time"] & merged["same_storm"]]
    if len(timed):
        print("   lag in HOURS, the %d pairs where BOTH records carry a clock "
              "time:" % len(timed))
        print("      median %+.1f   mean %+.1f   IQR %+.1f to %+.1f   "
              "range %+.1f to %+.1f"
              % (timed["lag_hours"].median(), timed["lag_hours"].mean(),
                 timed["lag_hours"].quantile(0.25),
                 timed["lag_hours"].quantile(0.75),
                 timed["lag_hours"].min(), timed["lag_hours"].max()))
    else:
        print("   no pair has a clock time on both sides -- the peak file "
              "publishes")
        print("   dates only for most of this record, so Part 2 carries the "
              "sub-daily answer.")
    off = merged[~merged["same_storm"]]
    if len(off):
        print("   %d year(s) more than %d days apart, treated as different "
              "storms: %s" % (len(off), PEAK_MAX_SEPARATION_DAYS,
                              ", ".join("WY%d" % w for w in off["WY"])))


# ----------------------------------------------------------------------------
# Part 2 -- sub-daily
# ----------------------------------------------------------------------------

def find_events(series, n, min_separation_days, floor):
    """Largest independent peaks: take the max, blank its neighbourhood, repeat."""
    working = series.copy()
    events = []
    span = pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span)
                    & (working.index <= when + span)] = np.nan
    return sorted(events)


def xcorr_lag(a, b, max_lag_hours, step_hours):
    """Lag in hours that best aligns b to a, by correlation of the anomalies.

    Both series are de-meaned first, so the match is driven by the shape of the
    rise and fall rather than by the difference in size between a large river
    and a small one.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    good = np.isfinite(a) & np.isfinite(b)
    if good.sum() < 12:
        return np.nan, np.nan
    a = a - np.nanmean(a[good])
    b = b - np.nanmean(b[good])
    a[~np.isfinite(a)] = 0.0
    b[~np.isfinite(b)] = 0.0
    max_shift = int(round(max_lag_hours / step_hours))
    best, best_r = np.nan, -np.inf
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            x, y = a[shift:], b[:len(b) - shift] if shift else b
        else:
            x, y = a[:len(a) + shift], b[-shift:]
        m = min(len(x), len(y))
        if m < 12:
            continue
        x, y = x[:m], y[:m]
        denom = np.sqrt(np.sum(x * x) * np.sum(y * y))
        if denom <= 0:
            continue
        r = float(np.sum(x * y) / denom)
        if r > best_r:
            best_r, best = r, shift
    if not np.isfinite(best_r):
        return np.nan, np.nan
    # NEGATED, and this is not cosmetic. The loop pairs a[shift + i] with b[i],
    # so if b genuinely arrives L hours AFTER a the best correlation is found
    # at shift = -L. Returning shift directly reported the Coweeman as leading
    # when it lags, reversing the conclusion the whole script exists to reach.
    # Verified against synthetic pairs built with a known offset.
    return -best * step_hours, best_r


def compare_subdaily(cow, cas):
    """Per-event lag between the two sub-daily records."""
    if cow is None or cas is None:
        return None
    first = max(cow.index.min(), cas.index.min())
    last = min(cow.index.max(), cas.index.max())
    if pd.isna(first) or pd.isna(last) or first >= last:
        print("   no overlap between the two sub-daily records")
        return None
    grid = pd.date_range(first.ceil(RESAMPLE), last.floor(RESAMPLE),
                         freq=RESAMPLE)
    cow_r = cow.resample(RESAMPLE).mean().reindex(grid)
    cas_r = cas.resample(RESAMPLE).mean().reindex(grid)
    step_hours = pd.Timedelta(RESAMPLE) / pd.Timedelta(hours=1)
    print("   overlap %s to %s (%d hours on a %s grid)"
          % (first.date(), last.date(), len(grid), RESAMPLE))

    events = find_events(cas_r, N_EVENTS, EVENT_MIN_SEPARATION_DAYS,
                         MIN_EVENT_CFS)
    print("   %d independent Castle Rock events above %s cfs"
          % (len(events), format(int(MIN_EVENT_CFS), ",")))
    half = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
    rows = []
    for when, value in events:
        window = slice(when - half, when + half)
        a, b = cas_r.loc[window], cow_r.loc[window]
        need = int(MIN_WINDOW_COVERAGE * len(a))
        if len(a) == 0 or a.notna().sum() < need or b.notna().sum() < need:
            continue
        t_cas, t_cow = a.idxmax(), b.idxmax()
        lag_peak = (t_cow - t_cas).total_seconds() / 3600.0
        lag_x, r = xcorr_lag(a.values, b.values, XCORR_MAX_LAG_HOURS,
                             step_hours)
        rows.append({"event_time": when, "cas_peak_cfs": float(a.max()),
                     "cow_peak_cfs": float(b.max()),
                     "cas_peak_time": t_cas, "cow_peak_time": t_cow,
                     "lag_peak_hours": lag_peak,
                     "lag_xcorr_hours": lag_x, "xcorr_r": r,
                     "ratio_cow_over_cas": float(b.max()) / float(a.max())
                     if a.max() else np.nan})
    return pd.DataFrame(rows) if rows else None


def report_subdaily(table):
    print("\n" + "=" * 78)
    print("PART 2 -- SUB-DAILY RECORDS")
    print("=" * 78)
    if table is None or not len(table):
        print("   no usable events")
        return
    print("   %d events with enough coverage in both records" % len(table))
    for column, label in (("lag_peak_hours", "peak to peak"),
                          ("lag_xcorr_hours", "cross-correlation")):
        values = table[column].dropna()
        if not len(values):
            continue
        print("   lag by %-18s (Coweeman minus Castle Rock), hours:" % label)
        print("      median %+.1f   mean %+.1f   IQR %+.1f to %+.1f   "
              "range %+.1f to %+.1f"
              % (values.median(), values.mean(), values.quantile(0.25),
                 values.quantile(0.75), values.min(), values.max()))
    agree = (table["lag_peak_hours"] - table["lag_xcorr_hours"]).abs()
    print("   the two methods differ by a median of %.1f hours; %d event(s) "
          "differ by more than 12"
          % (agree.median(), int((agree > 12).sum())))
    print("   Coweeman peak as a fraction of Castle Rock: median %.3f"
          % table["ratio_cow_over_cas"].median())
    lead = int((table["lag_xcorr_hours"] < 0).sum())
    print("   Coweeman leads Castle Rock in %d of %d events"
          % (lead, int(table["lag_xcorr_hours"].notna().sum())))


# ----------------------------------------------------------------------------

def plot_all(merged, sub, cow_sub, cas_sub, stem):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))

    ax = axes[0][0]
    if merged is not None and len(merged):
        same = merged[merged["same_storm"]]
        ax.scatter(same["WY"], same["lag_days"], s=34, color=C_COW,
                   edgecolor="0.2", lw=0.5, label="same storm")
        off = merged[~merged["same_storm"]]
        if len(off):
            ax.scatter(off["WY"], off["lag_days"], s=34, facecolor="none",
                       edgecolor="0.5", label="different storm")
        ax.axhline(0, color="0.3", lw=1.0, ls="--")
        ax.legend(fontsize=8)
    ax.set_title("Annual peak date difference by water year", fontsize=10)
    ax.set_xlabel("Water year", fontsize=9)
    ax.set_ylabel("Coweeman minus Castle Rock (days)", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[0][1]
    if sub is not None and len(sub):
        ax.scatter(sub["lag_peak_hours"], sub["lag_xcorr_hours"], s=40,
                   color=C_CAS, edgecolor="0.2", lw=0.5)
        lim = [min(sub["lag_peak_hours"].min(), sub["lag_xcorr_hours"].min()) - 3,
               max(sub["lag_peak_hours"].max(), sub["lag_xcorr_hours"].max()) + 3]
        ax.plot(lim, lim, color="0.4", lw=1.0, ls="--")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
    ax.set_title("Sub-daily lag: two methods compared", fontsize=10)
    ax.set_xlabel("Peak to peak (hours)", fontsize=9)
    ax.set_ylabel("Cross-correlation (hours)", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1][0]
    if sub is not None and len(sub):
        values = sub["lag_xcorr_hours"].dropna()
        if len(values):
            ax.hist(values, bins=max(6, min(20, len(values))), color=C_COW,
                    edgecolor="white")
            ax.axvline(values.median(), color="k", lw=1.6,
                       label="median %+.1f h" % values.median())
            ax.axvline(0, color="0.4", lw=1.0, ls="--")
            ax.legend(fontsize=8)
    ax.set_title("Distribution of sub-daily lag (cross-correlation)",
                 fontsize=10)
    ax.set_xlabel("Coweeman minus Castle Rock (hours)", fontsize=9)
    ax.set_ylabel("Events", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1][1]
    if sub is not None and len(sub) and cow_sub is not None:
        biggest = sub.loc[sub["cas_peak_cfs"].idxmax()]
        half = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
        window = slice(biggest["event_time"] - half,
                       biggest["event_time"] + half)
        a, b = cas_sub.loc[window], cow_sub.loc[window]
        if len(a) and len(b):
            ax.plot(a.index, a.values, color=C_CAS, lw=1.8,
                    label="Castle Rock")
            twin = ax.twinx()
            twin.plot(b.index, b.values, color=C_COW, lw=1.8,
                      label="Coweeman")
            twin.set_ylabel("Coweeman (cfs)", color=C_COW, fontsize=9)
            ax.set_ylabel("Castle Rock (cfs)", color=C_CAS, fontsize=9)
            ax.set_title("Largest common event: %s  (lag %+.1f h)"
                         % (biggest["event_time"].date(),
                            biggest["lag_xcorr_hours"]), fontsize=10)
            for label in ax.get_xticklabels():
                label.set_rotation(20)
                label.set_horizontalalignment("right")
    ax.grid(alpha=0.3)

    fig.suptitle("Coweeman River (%s / Ecology %s) against Cowlitz at Castle "
                 "Rock (%s)\npeak timing from the annual peak records and from "
                 "the sub-daily records"
                 % (COWEEMAN_USGS_SITE, ECOLOGY_SITE, CASTLE_ROCK_SITE),
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("%s.png" % stem, dpi=150)
    plt.close(fig)


def main():
    ensure_dirs()
    print("=" * 78)
    print("Coweeman vs Castle Rock peak timing")
    print("   cache : %s" % os.path.abspath(CACHE_DIR))
    print("=" * 78)

    print("\nDOWNLOAD")
    cow_peaks = fetch_usgs_peaks(COWEEMAN_USGS_SITE)
    cas_peaks = fetch_usgs_peaks(CASTLE_ROCK_SITE)
    cow_sub = fetch_ecology(ECOLOGY_SITE, ECOLOGY_FIRST_YEAR,
                            ECOLOGY_LAST_YEAR)
    cas_sub = None
    if cow_sub is not None:
        cas_sub = fetch_usgs_iv(CASTLE_ROCK_SITE,
                                str(cow_sub.index.min().date()),
                                str(cow_sub.index.max().date()))

    merged = compare_annual_peaks(cow_peaks, cas_peaks)
    report_annual(merged)
    if merged is not None:
        path = os.path.join(OUT_DIR, "coweeman_peak_timing.csv")
        merged.to_csv(path, index=False)
        print("   written: %s" % path)

    print("\nSUB-DAILY")
    sub = compare_subdaily(cow_sub, cas_sub) if cow_sub is not None else None
    report_subdaily(sub)
    if sub is not None:
        path = os.path.join(OUT_DIR, "coweeman_event_timing.csv")
        sub.to_csv(path, index=False)
        print("   written: %s" % path)

    if merged is not None or sub is not None:
        cow_r = (cow_sub.resample(RESAMPLE).mean()
                 if cow_sub is not None else None)
        cas_r = (cas_sub.resample(RESAMPLE).mean()
                 if cas_sub is not None else None)
        plot_all(merged, sub, cow_r, cas_r, PLOT_STEM)
        print("\nPlot : %s.png" % PLOT_STEM)

    print("\n" + "=" * 78)
    print("READ THE TWO PARTS TOGETHER")
    print("   Part 1 covers the long record but mostly at day resolution.")
    print("   Part 2 resolves hours but only over the Ecology period, and")
    print("   pairs two different agencies' gages. Agreement between them is")
    print("   the check; a disagreement means one of the two is not measuring")
    print("   what it appears to.")
    print("=" * 78)


main()
