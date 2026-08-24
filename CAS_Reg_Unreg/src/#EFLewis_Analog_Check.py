#EFLewis_Analog_Check.py
# -*- coding: utf-8 -*-
"""
Does the East Fork Lewis show the same high-AEP / low-AEP change in
coincident ratio that the Coweeman does?

WHY THIS EXISTS
    The downstream-confluence combination (#Coincident_TieredScaling.py)
    uses a two-tier factor: 0.80 for AEP>1%, 0.50 for AEP<=1%. The 0.80
    is well supported -- it is both the 2009-era East Fork Lewis
    analog-basin figure AND, independently, the Coweeman's own median in
    its 20,000-40,000 cfs bin (0.809, n=51). The 0.50 rests on the
    Coweeman's >60,000 cfs bin alone: n=8, several readings carrying
    non-good Ecology quality codes, and one (Dec 2015) known bad. See
    Downstream_Confluence_Notes.md sections 4-5.

    That is thin evidence for the number that drives the flood tail. The
    question this script answers: does an INDEPENDENT basin, gaged by a
    DIFFERENT agency, show the same drop-off at large Castle Rock events?
      - If yes, the tail drop is a real hydrologic effect and 0.50 is
        defensible for Arkansas and Ostrander too.
      - If no, the Coweeman's tail is more likely an artifact of its own
        rating-curve capping, and 0.80 flat becomes the better default.

    East Fork Lewis (USGS 14222500) is the right analog for the same
    reason the 2009 study picked it: 125 sq mi, close to the Coweeman's
    119, and it is the best-correlated neighbor gage per Duren (2015).
    It is 25 miles south-southeast, in a different basin -- a proxy for
    "how a small basin times against Castle Rock", not a Cowlitz
    tributary.

WHAT IS MEASURED
    Exactly what #Coweeman_RegPeak_Timing.py measures for the Coweeman,
    so the two are directly comparable:
      ratio = EF Lewis flow at the moment of the Castle Rock peak
              / EF Lewis's own peak within the event window
    Computed at BOTH the unregulated and the regulated Castle Rock peak
    timing, binned by the same Castle Rock magnitude bins.

NETWORK
    Needs one USGS instantaneous-values download (site 14222500), cached
    to CACHE_DIR on first run exactly like #Coweeman_Timing.py does, so a
    re-run is offline. THIS DOWNLOAD IS BLOCKED FROM THE CLAUDE SANDBOX
    (all usgs.gov hosts answer 403 at the network-policy layer, verified
    2026-08-24), so this script was written but never executed against
    real data. Run it locally. If the fetch fails it says so and exits
    rather than producing a half-answer.

    IV COVERAGE IS THE THING TO WATCH. The USGS instantaneous service
    generally starts around 2007. The 2009-era study cited 17 events over
    WY1996-2007, so it had sub-daily data this service may not return.
    The script prints the period it actually got and the per-bin counts,
    and warns if the large-event bin is too thin to conclude anything --
    a small-n answer here is no better than the small-n Coweeman answer
    it is meant to check.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pydsstools.heclib.dss import HecDss

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
EFLEWIS_SITE = "14222500"        # East Fork Lewis River near Heisson, WA
IV_PARAMETER = "00060"           # discharge, cfs
IV_START = "1990-01-01"          # ask wide; the service returns what it has
IV_END = "2020-12-31"

CACHE_DIR = r"../data/eflewis"
DSS_PATH = r"../output/ResSim_WCM_RC.dss"
UNREG_PATHNAME = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"
REG_PATHNAME = "//CastleRock_NWS/Flow//1Hour/ResSim_WCM_RC/"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "eflewis_analog_check.csv")
PLOT_PNG = os.path.join(OUT_DIR, "eflewis_analog_check.png")

# Same selection and binning as #Coweeman_RegPeak_Timing.py -- changing
# either of these breaks the comparability that is the whole point.
RESAMPLE = "1h"
N_EVENTS = 80
EVENT_MIN_SEPARATION_DAYS = 7
MIN_EVENT_CFS = 20000.0
EVENT_WINDOW_HOURS = 48
REG_SEARCH_WINDOW_HOURS = 72
MIN_WINDOW_COVERAGE = 0.80
MAGNITUDE_BINS = [20000, 40000, 60000, 200000]
BIN_LABELS = ["20-40k", "40-60k", ">60k"]
MIN_N_TO_CONCLUDE = 6            # below this a bin is reported, not believed

LOCAL_STANDARD_OFFSET_HOURS = -8
HTTP_HEADERS = {"User-Agent": "Mozilla/5.0 (compatible; USACE-hydrology-script)"}

# Coweeman's own numbers, from #Coweeman_RegPeak_Timing.py's run, for the
# side-by-side. Update if that script's inputs change.
COWEEMAN_REF = {
    "20-40k": {"n": 51, "unreg": 0.809, "reg": 0.809},
    "40-60k": {"n": 18, "unreg": 0.762, "reg": 0.762},
    ">60k":   {"n":  8, "unreg": 0.520, "reg": 0.494},
}

C_EFL = "#4c8c4a"
C_COW = "#b7410e"

# ----------------------------------------------------------------------------


def ensure_dirs():
    for path in (CACHE_DIR, OUT_DIR):
        if path and not os.path.isdir(path):
            os.makedirs(path)


def to_naive_local(index):
    index = pd.DatetimeIndex(index)
    if index.tz is not None:
        offset = pd.Timedelta(hours=LOCAL_STANDARD_OFFSET_HOURS)
        index = index.tz_convert("UTC") + offset
        index = index.tz_localize(None)
    return index


def parse_rdb(text):
    """USGS RDB -> DataFrame. Line 2 is the format spec, not data."""
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


def http_text(url, timeout=180):
    import requests
    try:
        response = requests.get(url, timeout=timeout, headers=HTTP_HEADERS)
        if response.status_code != 200:
            print("      HTTP %d  %s" % (response.status_code, url))
            return None
        return response.text
    except Exception as exc:
        print("      failed: %s" % exc)
        return None


def fetch_iv(site, start, end):
    """Instantaneous discharge, cached. Same idiom as #Coweeman_Timing.py:
    try dataretrieval first, fall back to the plain RDB service."""
    path = os.path.join(CACHE_DIR, "usgs_iv_%s_%s_%s.csv" % (site, start, end))
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
            column = next((c for c in data.columns
                           if c.startswith(IV_PARAMETER) and not c.endswith("_cd")), None)
            if column:
                series = pd.to_numeric(data[column], errors="coerce").dropna()
    except ImportError:
        print("   IV %s: dataretrieval not installed, using the RDB service" % site)
    except Exception as exc:
        print("   IV %s: dataretrieval failed (%s), using the RDB service" % (site, exc))

    if series is None:
        url = ("https://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=%s"
               "&startDT=%s&endDT=%s&parameterCd=%s" % (site, start, end, IV_PARAMETER))
        print("   GET %s" % url)
        table = parse_rdb(http_text(url))
        if table is None:
            return None
        value_col = next((c for c in table.columns
                          if c.endswith("_" + IV_PARAMETER)), None)
        if value_col is None:
            value_col = next((c for c in table.columns
                              if IV_PARAMETER in c and not c.endswith("_cd")), None)
        if value_col is None or "datetime" not in table.columns:
            print("   IV %s: no discharge column in %s" % (site, list(table.columns)))
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


def first_stamp(ts):
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


def read_dss_series(path, pathname):
    version = 6 if open(path, "rb").read(16)[12] == 6 else 7
    dss = HecDss.Open(path, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        values[np.array(ts.nodata, dtype=bool)] = np.nan
        values[values <= -900.0] = np.nan
        step = pd.Timedelta(seconds=int(ts.interval))
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).dropna().sort_index()


def find_events(series, n, min_separation_days, floor):
    working = series.copy()
    events, span = [], pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span) & (working.index <= when + span)] = np.nan
    return sorted(events)


def main():
    ensure_dirs()

    print("Reading series...")
    efl = fetch_iv(EFLEWIS_SITE, IV_START, IV_END)
    if efl is None or not len(efl):
        raise SystemExit(
            "\nEast Fork Lewis IV download FAILED -- no cached copy either.\n"
            "This script needs one USGS fetch to do anything; it does not\n"
            "guess. All usgs.gov hosts are blocked from the Claude sandbox\n"
            "(403 at CONNECT), so run this on a machine with USGS access.\n"
            "Once the CSV is cached under %s the script runs offline."
            % os.path.abspath(CACHE_DIR))

    unreg = read_dss_series(DSS_PATH, UNREG_PATHNAME)
    reg = read_dss_series(DSS_PATH, REG_PATHNAME)
    print("   Unreg CasRk: %d hours, %s to %s"
          % (len(unreg), unreg.index.min().date(), unreg.index.max().date()))

    first = max(efl.index.min(), unreg.index.min(), reg.index.min())
    last = min(efl.index.max(), unreg.index.max(), reg.index.max())
    grid = pd.date_range(first.ceil(RESAMPLE), last.floor(RESAMPLE), freq=RESAMPLE)
    efl_h = efl.resample(RESAMPLE).mean().reindex(grid)
    unreg_h = unreg.resample(RESAMPLE).mean().reindex(grid)
    reg_h = reg.resample(RESAMPLE).mean().reindex(grid)
    print("   overlap    : %s to %s (%d hours)" % (first.date(), last.date(), len(grid)))
    if (last - first) < pd.Timedelta(days=365 * 5):
        print("   WARNING: under 5 years of overlap -- too short to say anything "
              "about the large-event bin.")

    events = find_events(unreg_h, N_EVENTS, EVENT_MIN_SEPARATION_DAYS, MIN_EVENT_CFS)
    print("   %d independent unregulated events above %s cfs"
          % (len(events), format(int(MIN_EVENT_CFS), ",")))

    half_cow = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
    half_reg = pd.Timedelta(hours=REG_SEARCH_WINDOW_HOURS)
    rows = []
    for when, unreg_peak in events:
        window = slice(when - half_cow, when + half_cow)
        a, b = unreg_h.loc[window], efl_h.loc[window]
        need = int(MIN_WINDOW_COVERAGE * len(a))
        if len(a) == 0 or a.notna().sum() < need or b.notna().sum() < need:
            continue
        at_unreg = b.get(when, np.nan)
        if not np.isfinite(at_unreg):
            continue
        efl_peak = float(b.max())
        if efl_peak <= 0:
            continue

        c = reg_h.loc[slice(when - half_reg, when + half_reg)]
        if c.notna().sum() < int(MIN_WINDOW_COVERAGE * len(c)):
            continue
        reg_peak_time = c.idxmax()
        at_reg = efl_h.get(reg_peak_time, np.nan)
        if not np.isfinite(at_reg):
            continue

        rows.append({
            "event_time_unreg": when,
            "cas_unreg_peak_cfs": unreg_peak,
            "cas_reg_peak_cfs": float(c.max()),
            "reg_peak_time": reg_peak_time,
            "lag_reg_vs_unreg_hours": (reg_peak_time - when).total_seconds() / 3600.0,
            "efl_peak_cfs": efl_peak,
            "efl_at_unreg_peak_cfs": float(at_unreg),
            "efl_at_reg_peak_cfs": float(at_reg),
            "ratio_at_unreg_peak": float(at_unreg) / efl_peak,
            "ratio_at_reg_peak": float(at_reg) / efl_peak,
        })

    if not rows:
        raise SystemExit("No events survived the coverage filter -- check the "
                          "overlap period printed above.")

    table = pd.DataFrame(rows)
    table["magnitude_bin"] = pd.cut(table["cas_unreg_peak_cfs"], bins=MAGNITUDE_BINS,
                                     labels=BIN_LABELS, right=False)
    table.to_csv(OUT_CSV, index=False)
    print("Wrote %s (%d events)" % (OUT_CSV, len(table)))

    print()
    print("EAST FORK LEWIS vs COWEEMAN -- coincident ratio by Castle Rock magnitude")
    print("%-8s %-26s  %-26s" % ("", "East Fork Lewis (this run)", "Coweeman (reference)"))
    print("%-8s %5s %8s %8s   %5s %8s %8s" % ("bin", "n", "unreg", "reg", "n", "unreg", "reg"))
    verdict_rows = []
    for label in BIN_LABELS:
        g = table[table["magnitude_bin"] == label]
        ref = COWEEMAN_REF.get(label, {})
        if len(g) == 0:
            print("%-8s %5d %8s %8s   %5s %8s %8s"
                  % (label, 0, "--", "--", ref.get("n", "--"),
                     ref.get("unreg", "--"), ref.get("reg", "--")))
            continue
        eu, er = g["ratio_at_unreg_peak"].median(), g["ratio_at_reg_peak"].median()
        flag = "" if len(g) >= MIN_N_TO_CONCLUDE else "  <- too few to conclude"
        print("%-8s %5d %8.3f %8.3f   %5s %8s %8s%s"
              % (label, len(g), eu, er, ref.get("n", "--"),
                 ref.get("unreg", "--"), ref.get("reg", "--"), flag))
        verdict_rows.append((label, len(g), eu, er))

    print()
    tail = [r for r in verdict_rows if r[0] == ">60k"]
    common = [r for r in verdict_rows if r[0] == "20-40k"]
    if tail and common and tail[0][1] >= MIN_N_TO_CONCLUDE:
        drop = 100 * (1 - tail[0][2] / common[0][2])
        print("READ: East Fork Lewis drops %.0f%% from its 20-40k bin to its >60k bin "
              "(%.3f -> %.3f)." % (drop, common[0][2], tail[0][2]))
        print("      The Coweeman drops %.0f%% over the same bins (%.3f -> %.3f)."
              % (100 * (1 - COWEEMAN_REF[">60k"]["unreg"] / COWEEMAN_REF["20-40k"]["unreg"]),
                 COWEEMAN_REF["20-40k"]["unreg"], COWEEMAN_REF[">60k"]["unreg"]))
        print("      A similar drop in BOTH basins supports 0.50 as a real hydrologic")
        print("      effect. A flat East Fork Lewis points at Coweeman rating-curve")
        print("      capping instead, and argues for 0.80 throughout.")
    else:
        print("READ: not enough large events in the East Fork Lewis overlap to settle")
        print("      the question. This is the same small-n problem the check was")
        print("      meant to escape -- report it as unresolved, do not split the")
        print("      difference.")

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.scatter(table["cas_unreg_peak_cfs"], table["ratio_at_unreg_peak"],
               color=C_EFL, s=32, alpha=0.75, label="East Fork Lewis (this run)")
    for label in BIN_LABELS:
        g = table[table["magnitude_bin"] == label]
        if len(g):
            ax.hlines(g["ratio_at_unreg_peak"].median(),
                      g["cas_unreg_peak_cfs"].min(), g["cas_unreg_peak_cfs"].max(),
                      color=C_EFL, lw=2.5)
        ref = COWEEMAN_REF.get(label)
        if ref:
            lo, hi = MAGNITUDE_BINS[BIN_LABELS.index(label)], MAGNITUDE_BINS[BIN_LABELS.index(label) + 1]
            ax.hlines(ref["unreg"], lo, hi, color=C_COW, lw=2.5, ls="--")
    ax.plot([], [], color=C_COW, lw=2.5, ls="--", label="Coweeman bin median (reference)")
    ax.axhline(0.80, color="gray", lw=1, ls=":", label="0.80 tier factor")
    ax.axhline(0.50, color="gray", lw=1, ls="-.", label="0.50 tier factor")
    ax.set_xscale("log")
    ax.set_xlabel("Unregulated Castle Rock peak (cfs)")
    ax.set_ylabel("Tributary flow at Castle Rock peak / its own peak")
    ax.set_title("East Fork Lewis analog check: does the tail drop-off appear\n"
                  "in an independent basin, or only in the Coweeman?")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


if __name__ == "__main__":
    main()
