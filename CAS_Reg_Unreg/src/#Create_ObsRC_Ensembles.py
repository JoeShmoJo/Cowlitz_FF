#Create_ObsRC_Ensembles.py
# -*- coding: utf-8 -*-
"""
Build ResSim ensemble members centred on each water year's Castle Rock peak.

One member per water year. Each member starts at the BASE OF THE RISING LIMB of
the inflow hydrograph that produced that year's peak at Castle Rock, and runs
WINDOW_DAYS (default one month) -- long enough to contain the flood.

Three records are written per member:

    //MOSSYROCK/FLOW-IN//1HOUR/C:00000N|/     main window
    //CASTLE ROCK/FLOW-LOCAL//1HOUR/C:00000N|/ main window
    //MOS/ELEV//1HOUR/C:00000N|/               main window + LOOKBACK_DAYS earlier

The elevation is the OBSERVED pool from //MOS/ELEV//1DAY/USGS/, resampled to
hourly. It is a lookback record, so it starts LOOKBACK_DAYS before the flows.

Peak timing comes from the REGULATED Castle Rock record produced by the WCM rule
curve ResSim run (//CASTLEROCK_NWS/FLOW//1HOUR/ResSim_WCM_RC/), reassembled to a
period-of-record series by #Extract_Ensemble_To_Timeseries.py. The rising limb is
then found on the INFLOW hydrograph feeding that peak, so reservoir attenuation
between inflow and regulated outflow is handled by looking back from the
regulated peak rather than assuming the two coincide.

That regulated record only spans 01 Oct -> 01 May each year, so annual peaks are
by construction flood-season peaks. Water years with no regulated data are
skipped and listed.

Base detection smooths ONLY for detection. The values written out are the raw
volume-corrected inflows, bounces intact.

If MANUAL_STARTS_CSV exists its start dates are used verbatim, which is the
escape hatch when the automatic pick is wrong. The automatic picks are always
written to AUTO_STARTS_CSV in the same format, so it can be copied, edited and
fed back in.

Outputs: the ensemble DSS, a mapping CSV for reassembling results, the auto/edit
starts CSV, a per-member diagnostic CSV, and paged window plots.
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
from datetime import datetime, timedelta
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
IN_DSS = r"../output/ResSimInflows.dss"
IN_DSS_VERSION = 6
OBS_DSS = r"../../CAS_Unreg_FF/data/obsData.dss"
OBS_DSS_VERSION = 6

# EXTERNAL: requires the ResSim watershed (not in this repository)
OUT_DSS = r"C:\Projects\Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis\watershed\NWP_CowlitzLewis\shared\ensemble_obs_rc.dss"
OUT_DSS_VERSION = 6

MAPPING_CSV = r"../output/ensemble_obs_rc_mapping.csv"
AUTO_STARTS_CSV = r"../output/ensemble_obs_rc_starts_auto.csv"
MANUAL_STARTS_CSV = r"../data/ensemble_obs_rc_starts_manual.csv"   # used if present
DIAG_CSV = r"../output/diagnostics/ensemble_obs_rc_events.csv"
PLOT_STEM = r"../output/diagnostics/ensemble_obs_rc_windows"

PATH_MOS_IN = "//MOSSYROCK/FLOW-IN/*/1HOUR/FOR_RESSIM/"
PATH_CAS_LOCAL = "//CASTLE ROCK/FLOW-LOCAL/*/1HOUR/FOR_RESSIM/"
PATH_MOS_ELEV_DAILY = "//MOS/ELEV/*/1DAY/USGS/"

# REGULATED Castle Rock flow -- annual peaks are taken from this record.
# Produced by running #Extract_Ensemble_To_Timeseries.py on the WCM_RC simulation.
PEAK_DSS = r"../output/ResSim_WCM_RC.dss"
PEAK_DSS_VERSION = 6
PEAK_PATH = "//CASTLEROCK_NWS/FLOW/*/1HOUR/ResSim_WCM_RC/"
# Fail loudly rather than silently timing off the wrong hydrograph
FALLBACK_TO_INFLOW_SUM = False

# Hydrograph the rising limb is measured on: "SUM" (MOS inflow + CAS local)
# or "MOS" (Mossyrock inflow only)
LIMB_SERIES = "SUM"

WINDOW_DAYS = 31            # member length
LOOKBACK_DAYS = 1           # elevation record starts this much earlier
WATER_YEAR_START_MONTH = 10

# --- base-of-rising-limb detection ---
SMOOTH_HOURS = 12           # centred smoothing, FOR DETECTION ONLY
MAX_LOOKBACK_DAYS = 20      # how far back of the peak to hunt for the base
BASE_TOL = 0.05             # base = last time before peak within this fraction
                            # of the rise above the trough (0.05 = 5%)
MIN_DAYS_AFTER_PEAK = 7     # peak must sit at least this far from the window end
MIN_HOURS_BEFORE_PEAK = 6   # and at least this far from the window start

FIRST_WATER_YEAR = None     # None = every year with full coverage
LAST_WATER_YEAR = None
MAX_MISSING_HOURS = 24         # skip a year only if gaps exceed this; the rest
                               # are carry-forward filled

# Synthetic calendar the members are labelled on (any non-leap year)
ENS_LABEL_START = datetime(1999, 10, 1, 0, 0)   # hour-beginning

ELEV_DAILY_TO_HOURLY = "interpolate"   # "interpolate" or "step"

# How much observed pool to write.
#   "lookback" : LOOKBACK_DAYS ending AT the simulation start -- an INITIAL
#                condition, so ResSim's rules control the pool from there on.
#   "full"     : the lookback plus the whole window. Only use this if you want
#                the observed pool imposed across the simulation, which defeats
#                the purpose of letting the release rules operate.
ELEV_EXTENT = "lookback"
CLIP_NEGATIVE_FLOW = True

SENTINEL = -901.0
PLOTS_PER_PAGE = 24

# ----------------------------------------------------------------------------


def first_stamp(ts):
    """First timestamp of a DSS series, across pydsstools versions."""
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M", "%d%b%Y %H%M%S", "%d%b%Y %H%M",
                "%d %B %Y %H:%M:%S", "%d %B %Y %H:%M"):
        try:
            return pd.Timestamp(datetime.strptime(text, fmt))
        except ValueError:
            continue
    return pd.Timestamp(text)


def series_step(ts, pathname):
    """Time step of a DSS regular series. ts.interval is in seconds."""
    seconds = int(getattr(ts, "interval", 0) or 0)
    if seconds > 0:
        return pd.Timedelta(seconds=seconds)
    e_part = pathname.split("/")[5].upper()
    lookup = {"1MIN": "1min", "15MIN": "15min", "30MIN": "30min",
              "1HOUR": "1h", "6HOUR": "6h", "12HOUR": "12h", "1DAY": "1D"}
    return pd.Timedelta(lookup.get(e_part, "1h"))


def read_dss_series(dss_file, pathname, version):
    """Read a DSS regular time series into a Series on period-BEGINNING labels."""
    dss = HecDss.Open(dss_file, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        nodata = np.array(ts.nodata, dtype=bool)
        values[nodata] = np.nan
        values[values <= -900.0] = np.nan
        step = series_step(ts, pathname)
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index()


def daily_to_hourly(daily, how):
    """Daily elevation to hourly. Daily values are labelled on their own day."""
    hourly_index = pd.date_range(daily.index[0],
                                 daily.index[-1] + pd.Timedelta(hours=23), freq="h")
    if how == "step":
        return daily.reindex(hourly_index, method="ffill")
    # interpolate between midday-to-midday so a daily mean sits mid-day
    anchored = daily.copy()
    anchored.index = anchored.index + pd.Timedelta(hours=12)
    return anchored.reindex(anchored.index.union(hourly_index)) \
                   .interpolate(method="time") \
                   .reindex(hourly_index)


def water_year(stamp):
    return stamp.year + (1 if stamp.month >= WATER_YEAR_START_MONTH else 0)


def find_peak_and_base(peak_series, limb, smoothed, wy):
    """Regulated annual peak at Castle Rock, and the inflow rising limb feeding it.

    peak_series is the REGULATED Castle Rock flow -- it sets the timing.
    limb / smoothed are the INFLOW hydrograph -- they set the base.
    """
    a = pd.Timestamp(wy - 1, WATER_YEAR_START_MONTH, 1)
    b = pd.Timestamp(wy, WATER_YEAR_START_MONTH, 1) - pd.Timedelta(hours=1)
    year = peak_series.loc[a:b].dropna()
    if year.empty:
        return None
    peak_time = year.idxmax()
    peak_value = float(year.max())

    look_from = peak_time - pd.Timedelta(days=MAX_LOOKBACK_DAYS)
    seg = smoothed.loc[look_from:peak_time].dropna()
    if len(seg) < 3:
        return None

    trough_value = float(seg.min())
    smoothed_peak = float(seg.max())
    threshold = trough_value + BASE_TOL * max(smoothed_peak - trough_value, 1e-9)
    # LAST time at or below the threshold before the peak -- the foot of the limb
    at_base = seg[seg <= threshold]
    base_time = at_base.index[-1] if len(at_base) else seg.idxmin()

    # keep the peak comfortably inside the window
    latest = peak_time - pd.Timedelta(hours=MIN_HOURS_BEFORE_PEAK)
    earliest = peak_time - pd.Timedelta(days=WINDOW_DAYS - MIN_DAYS_AFTER_PEAK)
    clamped = "none"
    if base_time > latest:
        base_time, clamped = latest, "late"
    if base_time < earliest:
        base_time, clamped = earliest, "early"

    inflow_window = limb.loc[base_time:peak_time + pd.Timedelta(days=3)].dropna()
    inflow_peak_time = inflow_window.idxmax() if len(inflow_window) else pd.NaT
    inflow_peak_cfs = float(inflow_window.max()) if len(inflow_window) else np.nan
    attenuation_hours = (int((peak_time - inflow_peak_time).total_seconds() // 3600)
                         if inflow_window.size else np.nan)

    return {"water_year": wy, "peak_time": peak_time, "peak_cfs": peak_value,
            "inflow_peak_time": inflow_peak_time,
            "inflow_peak_cfs": inflow_peak_cfs,
            "reg_peak_lag_hours": attenuation_hours,
            "base_time": base_time.floor("h"),
            "base_cfs": float(limb.get(base_time.floor("h"), np.nan)),
            "trough_smoothed_cfs": trough_value,
            "base_threshold_cfs": threshold,
            "lead_hours": int((peak_time - base_time).total_seconds() // 3600),
            "clamped": clamped}


def load_manual_starts(path):
    """Optional manual override: columns water_year, base_time."""
    if not path or not os.path.isfile(path):
        return None
    table = pd.read_csv(path)
    table.columns = [c.strip().lower() for c in table.columns]
    table["base_time"] = pd.to_datetime(table["base_time"])
    return dict(zip(table["water_year"].astype(int), table["base_time"]))


def build_container(pathname, values, start_time, units, data_type, interval_min):
    """TimeSeriesContainer built the way the installed pydsstools accepts."""
    try:
        tsc = TimeSeriesContainer(pathname, len(values), interval_min,
                                  values=list(values), start_time=start_time,
                                  data_units=units, data_type=data_type)
    except TypeError:
        tsc = TimeSeriesContainer()
        tsc.pathname = pathname
        tsc.startDateTime = start_time
        tsc.numberValues = len(values)
        tsc.interval = interval_min
        tsc.values = list(values)
        tsc.units = units
        tsc.type = data_type
    return tsc


def slice_window(series, start, n_hours):
    """Values over the window; returns array, missing count, and the index."""
    idx = pd.date_range(start, periods=n_hours, freq="h")
    vals = series.reindex(idx)
    n_missing = int(vals.isna().sum())
    return vals.values, n_missing, idx


def to_dss_values(values):
    filled = pd.Series(values).ffill()
    return np.where(np.isfinite(filled.values), filled.values, SENTINEL)


def fmt_dss(dt):
    return pd.Timestamp(dt).strftime("%d%b%Y %H%M").upper()


def d_part(start, n_hours):
    return "%s - %s" % (fmt_dss(start),
                        fmt_dss(pd.Timestamp(start) + pd.Timedelta(hours=n_hours - 1)))


def plot_windows(events, peak_series, limb, mos, cas, elev_hourly, stem):
    """Paged panels: hydrograph, detected base, peak, and the member window."""
    pages = int(np.ceil(len(events) / float(PLOTS_PER_PAGE)))
    ncol = 4
    for page in range(pages):
        chunk = events[page * PLOTS_PER_PAGE:(page + 1) * PLOTS_PER_PAGE]
        nrow = int(np.ceil(len(chunk) / float(ncol)))
        fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 3.2 * nrow),
                                 squeeze=False)
        for k in range(nrow * ncol):
            ax = axes[k // ncol][k % ncol]
            if k >= len(chunk):
                ax.axis("off")
                continue
            ev = chunk[k]
            base = pd.Timestamp(ev["base_time"])
            end = base + pd.Timedelta(days=WINDOW_DAYS) - pd.Timedelta(hours=1)
            view_a = base - pd.Timedelta(days=10)
            view_b = end + pd.Timedelta(days=3)

            ax.plot(limb.loc[view_a:view_b].index, limb.loc[view_a:view_b].values,
                    color="0.25", lw=1.2)
            ax.plot(mos.loc[view_a:view_b].index, mos.loc[view_a:view_b].values,
                    color="#2c7fb8", lw=0.8)
            ax.plot(cas.loc[view_a:view_b].index, cas.loc[view_a:view_b].values,
                    color="#4c9a2a", lw=0.8)
            reg = peak_series.loc[view_a:view_b]
            if len(reg):
                ax.plot(reg.index, reg.values, color="#e67e22", lw=1.5)
            ax.axvspan(base, end, color="#c0392b", alpha=0.10)
            ax.axvline(base, color="#c0392b", lw=1.4)
            ax.axvline(pd.Timestamp(ev["peak_time"]), color="#e67e22", lw=1.2, ls="--")
            if pd.notna(ev.get("inflow_peak_time")):
                ax.axvline(pd.Timestamp(ev["inflow_peak_time"]), color="0.45",
                           lw=1.0, ls=":")
            ax.axhline(ev["base_threshold_cfs"], color="0.6", lw=0.7, ls=":")
            ax.set_title("WY%d  peak %.0f cfs  lead %.1f d%s"
                         % (ev["water_year"], ev["peak_cfs"], ev["lead_hours"] / 24.0,
                            "" if ev["clamped"] == "none" else
                            "  CLAMPED %s" % ev["clamped"]),
                         fontsize=9,
                         color="k" if ev["clamped"] == "none" else "#c0392b")
            ax.tick_params(labelsize=6)
            for label in ax.get_xticklabels():
                label.set_rotation(20)
                label.set_horizontalalignment("right")
            ax.set_ylabel("Flow (cfs)", fontsize=7)
            ax.grid(alpha=0.25)

            ax2 = ax.twinx()
            ev_elev = elev_hourly.loc[view_a:view_b]
            ax2.plot(ev_elev.index, ev_elev.values, color="#8e44ad", lw=1.0, ls="-.")
            ax2.set_ylabel("Pool (ft)", fontsize=7)
            ax2.tick_params(labelsize=6)

        handles = [
            Line2D([], [], color="0.25", lw=1.2, label="Inflow hydrograph (MOS in + local)"),
            Line2D([], [], color="#2c7fb8", lw=0.8, label="Mossyrock inflow"),
            Line2D([], [], color="#4c9a2a", lw=0.8, label="Castle Rock local"),
            Line2D([], [], color="#e67e22", lw=1.5, label="REGULATED Castle Rock (WCM_RC)"),
            Line2D([], [], color="#c0392b", lw=1.4, label="Detected base / window start"),
            Line2D([], [], color="#e67e22", lw=1.2, ls="--", label="Regulated annual peak"),
            Line2D([], [], color="0.45", lw=1.0, ls=":", label="Inflow peak"),
            Line2D([], [], color="0.6", lw=0.7, ls=":", label="Base threshold"),
            Line2D([], [], color="#8e44ad", lw=1.0, ls="-.", label="Observed pool elevation"),
        ]
        fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8.5,
                   frameon=False)
        fig.suptitle("Ensemble windows, page %d of %d -- base of the rising limb "
                     "feeding each water year's REGULATED Castle Rock peak"
                     % (page + 1, pages),
                     fontsize=12)
        fig.tight_layout(rect=[0, 0.045, 1, 0.97])
        fig.savefig("%s_page%d.png" % (stem, page + 1), dpi=140)
        plt.close(fig)


def plot_overview(events, stem):
    """Lead time and peak magnitude across years, to spot odd picks fast."""
    frame = pd.DataFrame(events)
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    x = frame["water_year"].values
    colors = ["#c0392b" if c != "none" else "#2c7fb8" for c in frame["clamped"]]
    axes[0].bar(x, frame["peak_cfs"], color="0.7")
    axes[0].set_ylabel("Regulated annual peak (cfs)")
    axes[0].set_title("Regulated Castle Rock annual peak per water year "
                      "(red = window clamped to contain the peak)")
    axes[0].grid(axis="y", alpha=0.3)
    axes[1].bar(x, frame["lead_hours"] / 24.0, color=colors)
    axes[1].axhline(WINDOW_DAYS - MIN_DAYS_AFTER_PEAK, color="k", lw=0.9, ls="--")
    axes[1].set_ylabel("Base to peak (days)")
    axes[1].set_xlabel("Water year")
    axes[1].grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig("%s_overview.png" % stem, dpi=140)
    plt.close(fig)


def main():
    for path in (os.path.dirname(MAPPING_CSV), os.path.dirname(DIAG_CSV),
                 os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    mos = read_dss_series(IN_DSS, PATH_MOS_IN, IN_DSS_VERSION)
    cas = read_dss_series(IN_DSS, PATH_CAS_LOCAL, IN_DSS_VERSION)
    if CLIP_NEGATIVE_FLOW:
        n_neg = int((mos < 0).sum() + (cas < 0).sum())
        mos = mos.clip(lower=0.0)
        cas = cas.clip(lower=0.0)
    else:
        n_neg = 0

    elev_daily = read_dss_series(OBS_DSS, PATH_MOS_ELEV_DAILY, OBS_DSS_VERSION)
    elev_hourly = daily_to_hourly(elev_daily.dropna(), ELEV_DAILY_TO_HOURLY)

    joint = mos.index.union(cas.index)
    inflow_sum = (mos.reindex(joint).fillna(0.0) + cas.reindex(joint).fillna(0.0))
    inflow_sum = inflow_sum[mos.notna().reindex(joint, fill_value=False)
                            | cas.notna().reindex(joint, fill_value=False)]
    limb = inflow_sum if LIMB_SERIES == "SUM" else mos.dropna()
    smoothed = limb.rolling(SMOOTH_HOURS, center=True, min_periods=1).mean()

    # REGULATED Castle Rock flow sets the peak timing
    try:
        peak_series = read_dss_series(PEAK_DSS, PEAK_PATH, PEAK_DSS_VERSION).dropna()
        peak_source = "%s  %s" % (PEAK_DSS, PEAK_PATH)
    except Exception as exc:
        if not FALLBACK_TO_INFLOW_SUM:
            print("ERROR: could not read the regulated Castle Rock record.")
            print("   %s" % exc)
            print("   %s :: %s" % (PEAK_DSS, PEAK_PATH))
            print("   Run #Extract_Ensemble_To_Timeseries.py on the WCM_RC simulation")
            print("   first -- that is what produces this record. Set")
            print("   FALLBACK_TO_INFLOW_SUM = True only if you deliberately want")
            print("   UNREGULATED peak timing instead.")
            return
        print("WARNING: regulated record unavailable, timing off the inflow sum")
        peak_series = inflow_sum
        peak_source = "FALLBACK: unregulated inflow sum"

    n_hours = WINDOW_DAYS * 24
    if ELEV_EXTENT == "lookback":
        elev_hours = LOOKBACK_DAYS * 24 + 1      # ends ON the simulation start hour
    else:
        elev_hours = n_hours + LOOKBACK_DAYS * 24
    ens_start = pd.Timestamp(ENS_LABEL_START)
    elev_ens_start = ens_start - pd.Timedelta(days=LOOKBACK_DAYS)
    flow_d = d_part(ens_start, n_hours)
    elev_d = d_part(elev_ens_start, elev_hours)
    flow_write_start = (ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()
    elev_write_start = (elev_ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()

    manual = load_manual_starts(MANUAL_STARTS_CSV)
    years = sorted(set(water_year(t) for t in peak_series.index))
    if FIRST_WATER_YEAR:
        years = [y for y in years if y >= FIRST_WATER_YEAR]
    if LAST_WATER_YEAR:
        years = [y for y in years if y <= LAST_WATER_YEAR]

    print("=" * 78)
    print("Window        : %d days (%d hours) from the base of the rising limb"
          % (WINDOW_DAYS, n_hours))
    print("Elevation     : observed daily pool -> hourly (%s), %d day lookback, "
          "extent=%s (%d hours)"
          % (ELEV_DAILY_TO_HOURLY, LOOKBACK_DAYS, ELEV_EXTENT, elev_hours))
    print("Flow D-part   : %s" % flow_d)
    print("Elev D-part   : %s" % elev_d)
    print("Peak timing   : %s" % peak_source)
    print("Rising limb   : %s" % ("MOS inflow + CAS local" if LIMB_SERIES == "SUM"
                                  else "Mossyrock inflow only"))
    print("Start dates   : %s" % ("MANUAL from %s (%d years)" % (MANUAL_STARTS_CSV, len(manual))
                                  if manual else "auto-detected"))
    if CLIP_NEGATIVE_FLOW and n_neg:
        print("Clipped %d negative flow values to zero" % n_neg)
    print("=" * 78)

    events, skipped = [], []
    for wy in years:
        ev = find_peak_and_base(peak_series, limb, smoothed, wy)
        if ev is None:
            skipped.append((wy, "no regulated Castle Rock data in this water year"))
            continue
        if manual and wy in manual:
            ev["base_time"] = pd.Timestamp(manual[wy]).floor("h")
            ev["lead_hours"] = int((pd.Timestamp(ev["peak_time"])
                                    - ev["base_time"]).total_seconds() // 3600)
            ev["clamped"] = "manual"
        base = pd.Timestamp(ev["base_time"])
        _, mos_miss, _ = slice_window(mos, base, n_hours)
        _, cas_miss, _ = slice_window(cas, base, n_hours)
        _, elev_miss, _ = slice_window(elev_hourly, base - pd.Timedelta(days=LOOKBACK_DAYS),
                                       elev_hours)
        ev["mos_missing_hours"] = mos_miss
        ev["cas_missing_hours"] = cas_miss
        ev["elev_missing_hours"] = elev_miss
        ev["peak_in_window"] = bool(base <= pd.Timestamp(ev["peak_time"])
                                    <= base + pd.Timedelta(hours=n_hours - 1))
        if max(mos_miss, cas_miss) > MAX_MISSING_HOURS:
            skipped.append((wy, "missing %d MOS / %d CAS hours" % (mos_miss, cas_miss)))
            continue
        events.append(ev)

    if not events:
        print("No members built.")
        return

    mapping_rows = []
    with HecDss.Open(OUT_DSS, version=OUT_DSS_VERSION) as dst:
        for member, ev in enumerate(events, start=1):
            base = pd.Timestamp(ev["base_time"])
            f_part = "C:%06d|" % member
            mos_v, _, _ = slice_window(mos, base, n_hours)
            cas_v, _, _ = slice_window(cas, base, n_hours)
            elev_v, _, _ = slice_window(elev_hourly,
                                        base - pd.Timedelta(days=LOOKBACK_DAYS),
                                        elev_hours)
            for parts, vals, units, dpart, wstart in [
                    (("", "MOSSYROCK", "FLOW-IN"), mos_v, "CFS", flow_d, flow_write_start),
                    (("", "CASTLE ROCK", "FLOW-LOCAL"), cas_v, "CFS", flow_d, flow_write_start),
                    (("", "MOS", "ELEV"), elev_v, "FEET", elev_d, elev_write_start)]:
                pathname = "/%s/%s/%s/%s/1HOUR/%s/" % (parts[0], parts[1], parts[2],
                                                       dpart, f_part)
                dst.put_ts(build_container(pathname, to_dss_values(vals), wstart,
                                           units, "INST-VAL", 60))
            mapping_rows.append({
                "member": member, "ensemble_f_part": f_part,
                "water_year": ev["water_year"],
                "real_start": base,
                "real_end": base + pd.Timedelta(hours=n_hours - 1),
                "ensemble_start": ens_start, "hours": n_hours,
                "elev_real_start": base - pd.Timedelta(days=LOOKBACK_DAYS),
                "elev_ensemble_start": elev_ens_start, "elev_hours": elev_hours,
                "peak_time": ev["peak_time"], "peak_cfs": round(ev["peak_cfs"], 1),
                "lead_hours": ev["lead_hours"], "clamped": ev["clamped"],
                "start_pool_ft": round(float(elev_v[-1]), 2)})

    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(MAPPING_CSV, index=False)
    diag = pd.DataFrame(events)
    diag.to_csv(DIAG_CSV, index=False)
    diag[["water_year", "base_time"]].to_csv(AUTO_STARTS_CSV, index=False)

    plot_windows(events, peak_series, limb, mos, cas, elev_hourly, PLOT_STEM)
    plot_overview(events, PLOT_STEM)

    clamped = diag[diag["clamped"] != "none"]
    print("Members written : %d   records: %d" % (len(events), len(events) * 3))
    print("Mapping CSV     : %s" % MAPPING_CSV)
    print("Editable starts : %s  (copy to %s to override)"
          % (AUTO_STARTS_CSV, MANUAL_STARTS_CSV))
    print("Diagnostics CSV : %s" % DIAG_CSV)
    print("Plots           : %s_page*.png and %s_overview.png" % (PLOT_STEM, PLOT_STEM))
    print("-" * 78)
    pool = mapping["start_pool_ft"]
    print("Starting pool at event onset: median %.1f ft, min %.1f ft, max %.1f ft"
          % (pool.median(), pool.min(), pool.max()))
    print("Lead time base->peak: median %.1f d, min %.1f d, max %.1f d"
          % (diag["lead_hours"].median() / 24.0, diag["lead_hours"].min() / 24.0,
             diag["lead_hours"].max() / 24.0))
    print("Windows clamped to contain the peak: %d" % len(clamped))
    for _, r in clamped.head(12).iterrows():
        print("   WY%d  %s  lead %.1f d" % (r["water_year"], r["clamped"],
                                            r["lead_hours"] / 24.0))
    if not diag["peak_in_window"].all():
        print("*** %d windows do NOT contain their peak -- inspect the plots ***"
              % int((~diag["peak_in_window"]).sum()))
    if skipped:
        print("Skipped %d water years:" % len(skipped))
        for wy, why in skipped[:12]:
            print("   WY%d  %s" % (wy, why))


main()
