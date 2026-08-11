#Critical_Duration_Adjusted.py
# -*- coding: utf-8 -*-
"""
Critical duration for Castle Rock, using the ADJUSTED regulated peak record.

QUESTION
--------
Which unregulated flow duration best explains the regulated peak at Castle
Rock? If the regulated peak tracks the 1-day unregulated flow the system is
peak-limited; if it tracks the 5-day the system is volume-limited and the
reservoir fills before the flood ends. The duration with the strongest fit is
the critical duration.

WHAT IS DIFFERENT FROM Critical_Duration_Correlation.py
-------------------------------------------------------
1. REGULATED PEAK is the ADJUSTED record from #Adjusted_Peak_Record.py, not
   the raw USGS peak. Every year is on a consistent rule-curve starting-pool
   basis, so year-to-year variation in the observed starting pool is removed
   from the scatter.

2. UNREGULATED DURATIONS come from the ResSim unregulated period-of-record
   run, read at 1HOUR resolution, rather than from the daily mass balance.
   That makes the true hourly peak available and puts all four durations on
   one consistent source. Despite "Ensemble--0" in its F-part, that record is
   a single POR simulation, not an ensemble.

3. DURATIONS ARE EVENT-BASED, not water-year maxima: each duration is the
   maximum N-hour mean within +/- EVENT_WINDOW_DAYS of the regulated peak, so
   every duration describes the storm that produced that peak.

4. FITS ARE LOG-LOG. reg = a * unreg^b, fitted as log10(reg) on log10(unreg).
   A power law holds shape across magnitudes, which a linear fit does not --
   and this relationship is wanted for LARGE events, so it must not be
   dominated by the small ones. Linear statistics are reported alongside for
   reference only.

CROSS-CHECK
-----------
The ResSim unregulated durations are compared against the ones adopted in the
unregulated flow frequency study (unreg_durations_massbalance.csv and the
hourly holdout in wy_peak_records.csv). These are independent estimates of the
same quantity -- ResSim routing versus a daily storage mass balance -- so
systematic disagreement is worth knowing about before either is relied on.

CAVEAT
------
Unlike the mass-balance record, the ResSim unregulated flow is not built from
the regulated record, so this pairing is not self-correlated by construction.
The adjusted regulated peak does contain a ResSim-derived term, which is a
weaker version of the same concern -- the adjustment, not the peak itself.

OUTPUTS (../output and ../output/diagnostics)
    critical_duration_adjusted_dataset.csv       one row per water year
    critical_duration_adjusted_fits.csv          fit statistics per duration
    critical_duration_adjusted_scatter.png       log-log scatter with fits
    critical_duration_adjusted_summary.png       r-squared vs duration
    unreg_source_comparison.csv / .png           ResSim vs mass-balance durations
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime
from pydsstools.heclib.dss import HecDss

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
# EXTERNAL: the ResSim unregulated period-of-record run (not in the repository).
# "Ensemble--0" in the F-part is a ResSim artifact; this is a single POR run.
UNREG_DSS = (r"C:\Projects\2026_Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis"
             r"\watershed\NWP_CowlitzLewis_ResSim4\rss\Unreg_POR_FIS\simulation.dss")
UNREG_PATH = "//CastleRock_NWS/Flow-UNREG/*/1Hour/Ensemble--0/"

ADJUSTED_PEAKS_CSV = r"../output/adjusted_peaks.csv"

# Cross-check sources, from the unregulated flow frequency study
MASSBAL_CSV = r"../../CAS_Unreg_FF/output/unreg_durations_massbalance.csv"
WY_PEAKS_CSV = r"../../CAS_Unreg_FF/output/wy_peak_records.csv"

OUT_DIR = r"../output"
DIAG_DIR = r"../output/diagnostics"

WATER_YEAR_START_MONTH = 10
FIRST_WY = 1974              # Mossyrock/Riffe regulation begins
EXCLUDE_WYS = []             # the adjusted record already screens WY1980

# Durations: (label, hours). "Peak" is the instantaneous hourly value.
DURATIONS = [("Peak (1-hr)", 1),
             ("1-Day", 24),
             ("3-Day", 72),
             ("5-Day", 120)]

# Each duration is the max N-hour mean within +/- this many days of the
# regulated peak, so every duration describes the same storm.
EVENT_WINDOW_DAYS = 5
# A duration needs at least this fraction of its window present to be used.
MIN_COVERAGE = 0.95

# Use only years whose regulated peak passed the adjusted-record screening.
REQUIRE_SCREEN_PASSED = True
# Restrict the fit to events at or above this regulated peak. None = all years.
# The relationship is wanted for large events; set this if the small years are
# distorting the fit.
LARGE_EVENT_MIN_CFS = None

SENTINEL_MIN = -900.0

# ----------------------------------------------------------------------------


def dss_version(path):
    """DSS file version from the header: byte 12 is 6 for v6, 0 for v7.

    ResSim writes v7 in some runs and v6 in others, and pydsstools needs the
    version passed explicitly on Linux, so detect rather than assume.
    """
    with open(path, "rb") as handle:
        head = handle.read(16)
    if len(head) < 16 or head[:4] != b"ZDSS":
        return None
    return 6 if head[12] == 6 else 7


def first_stamp(ts):
    """First timestamp of a DSS series, across pydsstools versions."""
    first = next(iter(ts.times))
    if hasattr(first, "datetime"):
        return pd.Timestamp(first.datetime())
    text = str(getattr(ts, "startDateTime", None) or first).strip()
    # DSS uses midnight-as-2400: "01Oct1973 24:00:00" means 02Oct1973 00:00
    roll_day = False
    if " 24:" in text or text.endswith(" 2400"):
        text = text.replace(" 24:", " 00:").replace(" 2400", " 0000")
        roll_day = True
    for fmt in ("%d%b%Y %H:%M:%S", "%d%b%Y %H:%M", "%d%b%Y %H%M%S", "%d%b%Y %H%M",
                "%d %B %Y %H:%M:%S", "%d %B %Y %H:%M"):
        try:
            stamp = pd.Timestamp(datetime.strptime(text, fmt))
            return stamp + pd.Timedelta(days=1) if roll_day else stamp
        except ValueError:
            continue
    stamp = pd.Timestamp(text)
    return stamp + pd.Timedelta(days=1) if roll_day else stamp


def series_step(ts, pathname):
    """Time step of a DSS regular series. ts.interval is in seconds."""
    seconds = int(getattr(ts, "interval", 0) or 0)
    if seconds > 0:
        return pd.Timedelta(seconds=seconds)
    e_part = pathname.split("/")[5].upper()
    lookup = {"1MIN": "1min", "15MIN": "15min", "30MIN": "30min",
              "1HOUR": "1h", "6HOUR": "6h", "12HOUR": "12h", "1DAY": "1D"}
    return pd.Timedelta(lookup.get(e_part, "1h"))


def read_dss_series(dss_file, pathname):
    """Read a DSS regular time series into a Series on period-BEGINNING labels."""
    version = dss_version(dss_file)
    dss = HecDss.Open(dss_file, version=version) if version else HecDss.Open(dss_file)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        nodata = np.array(ts.nodata, dtype=bool)
        values[nodata] = np.nan
        values[values <= SENTINEL_MIN] = np.nan
        step = series_step(ts, pathname)
        index = pd.date_range(first_stamp(ts) - step, periods=len(values), freq=step)
    finally:
        dss.close()
    return pd.Series(values, index=index).sort_index()


def water_year(stamp):
    return stamp.year + (1 if stamp.month >= WATER_YEAR_START_MONTH else 0)


def event_duration(series, when, hours, half_width_days, min_coverage):
    """Max N-hour mean within +/- half_width_days of a time, plus its centre.

    Returns (value, centre_time, coverage). The rolling mean is right-labelled,
    so the window it covers ends at that label; the centre is reported for
    plotting and for checking the duration sits on the same storm.
    """
    pad = pd.Timedelta(days=half_width_days) + pd.Timedelta(hours=hours)
    window = series.loc[when - pad:when + pad]
    if len(window) == 0:
        return np.nan, pd.NaT, 0.0
    full = window.reindex(pd.date_range(window.index[0], window.index[-1], freq="h"))
    rolled = full.rolling(hours, min_periods=int(np.ceil(hours * min_coverage))).mean()
    keep = (rolled.index >= when - pd.Timedelta(days=half_width_days)) & \
           (rolled.index <= when + pd.Timedelta(days=half_width_days) +
            pd.Timedelta(hours=hours))
    rolled = rolled[keep].dropna()
    if len(rolled) == 0:
        return np.nan, pd.NaT, 0.0
    end = rolled.idxmax()
    start = end - pd.Timedelta(hours=hours - 1)
    coverage = float(full.loc[start:end].notna().mean())
    centre = start + (end - start) / 2
    return float(rolled.max()), centre, coverage


def fit_pair(x, y):
    """Linear and log-log fits of y on x, with the statistics that matter."""
    out = {"n": int(len(x))}
    if len(x) < 3:
        return out
    lin = stats.linregress(x, y)
    resid = y - (lin.slope * x + lin.intercept)
    out.update({"pearson_r": lin.rvalue, "r2": lin.rvalue ** 2,
                "p_value": lin.pvalue, "slope": lin.slope,
                "intercept": lin.intercept,
                "se_estimate_cfs": float(np.std(resid, ddof=2))})
    good = (x > 0) & (y > 0)
    if good.sum() >= 3:
        lx, ly = np.log10(x[good]), np.log10(y[good])
        log = stats.linregress(lx, ly)
        log_resid = ly - (log.slope * lx + log.intercept)
        # back-transformed scatter, as a multiplicative factor
        factor = 10 ** float(np.std(log_resid, ddof=2))
        out.update({"log_r2": log.rvalue ** 2, "log_exponent_b": log.slope,
                    "log_coeff_a": 10 ** log.intercept,
                    "log_p_value": log.pvalue,
                    "log_se_dex": float(np.std(log_resid, ddof=2)),
                    "log_se_factor": factor,
                    "spearman_rho": float(stats.spearmanr(x[good], y[good]).statistic),
                    "n_log": int(good.sum())})
    return out


def plot_scatter(data, fits, stem):
    """Log-log scatter with the fitted power law, one panel per duration."""
    labels = [d[0] for d in DURATIONS]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10))
    for k, label in enumerate(labels):
        ax = axes[k // 2][k % 2]
        col = "unreg_%s" % label.replace(" ", "_").replace("(", "").replace(")", "")
        sub = data.dropna(subset=[col, "reg_peak"])
        if len(sub) == 0:
            ax.axis("off")
            continue
        ax.loglog(sub[col], sub["reg_peak"], ls="none", marker="o", ms=5,
                  color="#2c7fb8", mec="0.3", mew=0.5)
        row = fits[fits["duration"] == label]
        if len(row) and np.isfinite(row["log_r2"].iloc[0]):
            a = row["log_coeff_a"].iloc[0]
            b = row["log_exponent_b"].iloc[0]
            xs = np.logspace(np.log10(sub[col].min()), np.log10(sub[col].max()), 50)
            ax.loglog(xs, a * xs ** b, color="#c0392b", lw=1.6)
            ax.set_title("%s\nreg = %.3g x unreg^%.3f   log r2 = %.3f   n = %d"
                         % (label, a, b, row["log_r2"].iloc[0], row["n_log"].iloc[0]),
                         fontsize=10)
        lim = [min(sub[col].min(), sub["reg_peak"].min()) * 0.8,
               max(sub[col].max(), sub["reg_peak"].max()) * 1.2]
        ax.plot(lim, lim, color="0.6", lw=0.9, ls="--")
        ax.set_xlabel("Unregulated %s (cfs)" % label)
        ax.set_ylabel("Adjusted regulated peak (cfs)")
        ax.grid(which="both", alpha=0.25)
    fig.suptitle("Critical duration -- adjusted regulated peak vs unregulated "
                 "event durations (log-log; dashed = 1:1)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("%s_scatter.png" % stem, dpi=150)
    plt.close(fig)


def plot_summary(fits, stem):
    """Fit quality against duration -- the critical duration is the maximum."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(fits))
    ax.plot(x, fits["log_r2"], color="#c0392b", lw=1.8, marker="o", ms=7,
            label="log-log r-squared")
    ax.plot(x, fits["r2"], color="0.55", lw=1.2, marker="s", ms=5, ls="--",
            label="linear r-squared (reference)")
    best = int(np.nanargmax(fits["log_r2"].values))
    ax.axvline(best, color="#16a085", lw=1.2, ls=":")
    ax.annotate("critical duration: %s" % fits["duration"].iloc[best],
                (best, fits["log_r2"].iloc[best]), xytext=(8, -14),
                textcoords="offset points", fontsize=10, color="#16a085")
    ax.set_xticks(x)
    ax.set_xticklabels(fits["duration"])
    ax.set_ylabel("r-squared")
    ax.set_xlabel("Unregulated flow duration")
    ax.set_title("Which unregulated duration explains the adjusted regulated peak?")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("%s_summary.png" % stem, dpi=150)
    plt.close(fig)


def plot_source_comparison(compare, stem):
    """ResSim unregulated durations against the mass-balance record."""
    pairs = [("unreg_Peak_1-hr", "mb_peak_1hr", "Peak (1-hr)"),
             ("unreg_1-Day", "mb_One_day", "1-Day"),
             ("unreg_3-Day", "mb_Three_Day", "3-Day"),
             ("unreg_5-Day", "mb_Five_Day", "5-Day")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for k, (a_col, b_col, label) in enumerate(pairs):
        ax = axes[k // 2][k % 2]
        if a_col not in compare or b_col not in compare:
            ax.axis("off")
            continue
        sub = compare.dropna(subset=[a_col, b_col])
        if len(sub) == 0:
            ax.axis("off")
            continue
        ax.loglog(sub[b_col], sub[a_col], ls="none", marker="o", ms=5,
                  color="#8e44ad", mec="0.3", mew=0.5)
        lim = [min(sub[a_col].min(), sub[b_col].min()) * 0.8,
               max(sub[a_col].max(), sub[b_col].max()) * 1.2]
        ax.plot(lim, lim, color="k", lw=1.0, ls="--")
        ratio = (sub[a_col] / sub[b_col]).median()
        ax.set_title("%s   n = %d   median ResSim/mass-balance = %.3f"
                     % (label, len(sub), ratio), fontsize=10)
        ax.set_xlabel("Mass-balance unregulated (cfs)")
        ax.set_ylabel("ResSim unregulated (cfs)")
        ax.grid(which="both", alpha=0.25)
    fig.suptitle("Unregulated durations: ResSim POR run vs the flow frequency "
                 "study's mass balance (dashed = 1:1)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("%s_source_comparison.png" % stem, dpi=150)
    plt.close(fig)


def main():
    for path in (OUT_DIR, DIAG_DIR):
        if not os.path.isdir(path):
            os.makedirs(path)

    if not os.path.isfile(UNREG_DSS):
        print("ERROR: the ResSim unregulated POR run was not found at")
        print("   %s" % UNREG_DSS)
        print("This file lives outside the repository. Correct UNREG_DSS.")
        return

    unreg = read_dss_series(UNREG_DSS, UNREG_PATH).dropna()
    print("=" * 78)
    print("Unregulated POR : %s .. %s  (%d hourly values)"
          % (unreg.index[0].date(), unreg.index[-1].date(), len(unreg)))
    print("   detected DSS version %s" % dss_version(UNREG_DSS))

    peaks = pd.read_csv(ADJUSTED_PEAKS_CSV, parse_dates=["t_usgs"])
    if REQUIRE_SCREEN_PASSED:
        peaks = peaks[peaks["screen_passed"]]
    peaks = peaks[(peaks["WY"] >= FIRST_WY) & (~peaks["WY"].isin(EXCLUDE_WYS))]
    print("Adjusted peaks  : %d water years (%d..%d)"
          % (len(peaks), peaks["WY"].min(), peaks["WY"].max()))
    print("Durations taken within +/- %d days of the regulated peak"
          % EVENT_WINDOW_DAYS)
    print("=" * 78)

    rows = []
    for _, row in peaks.iterrows():
        when = pd.Timestamp(row["t_usgs"])
        entry = {"WY": int(row["WY"]), "reg_peak": float(row["adjusted_peak"]),
                 "usgs_peak": float(row["usgs"]),
                 "adjustment": float(row["adjusted_peak"] - row["usgs"]),
                 "reg_peak_time": when}
        for label, hours in DURATIONS:
            key = label.replace(" ", "_").replace("(", "").replace(")", "")
            value, centre, coverage = event_duration(
                unreg, when, hours, EVENT_WINDOW_DAYS, MIN_COVERAGE)
            entry["unreg_%s" % key] = value
            entry["unreg_%s_centre" % key] = centre
            entry["unreg_%s_cov" % key] = coverage
            entry["unreg_%s_lag_hrs" % key] = (
                (centre - when).total_seconds() / 3600.0 if pd.notna(centre) else np.nan)
        rows.append(entry)
    data = pd.DataFrame(rows)

    if LARGE_EVENT_MIN_CFS:
        before = len(data)
        data = data[data["reg_peak"] >= LARGE_EVENT_MIN_CFS]
        print("Large-event filter >= %.0f cfs: %d of %d years kept"
              % (LARGE_EVENT_MIN_CFS, len(data), before))

    fit_rows = []
    for label, hours in DURATIONS:
        key = label.replace(" ", "_").replace("(", "").replace(")", "")
        sub = data.dropna(subset=["unreg_%s" % key, "reg_peak"])
        stats_row = {"duration": label, "duration_hours": hours,
                     "first_wy": int(sub["WY"].min()) if len(sub) else np.nan,
                     "last_wy": int(sub["WY"].max()) if len(sub) else np.nan}
        stats_row.update(fit_pair(sub["unreg_%s" % key].values,
                                  sub["reg_peak"].values))
        stats_row["median_reg_over_unreg"] = float(
            (sub["reg_peak"] / sub["unreg_%s" % key]).median()) if len(sub) else np.nan
        fit_rows.append(stats_row)
    fits = pd.DataFrame(fit_rows)

    data.to_csv(os.path.join(DIAG_DIR, "critical_duration_adjusted_dataset.csv"),
                index=False, float_format="%.2f")
    fits.to_csv(os.path.join(OUT_DIR, "critical_duration_adjusted_fits.csv"),
                index=False, float_format="%.5f")

    stem = os.path.join(DIAG_DIR, "critical_duration_adjusted")
    plot_scatter(data, fits, stem)
    plot_summary(fits, stem)

    print("\nFIT RESULTS (log-log: reg = a * unreg^b)")
    show = fits[["duration", "n_log", "log_r2", "log_exponent_b", "log_coeff_a",
                 "log_se_factor", "r2", "median_reg_over_unreg"]].copy()
    print(show.round(4).to_string(index=False))
    best = fits.loc[fits["log_r2"].idxmax()]
    print("\nCRITICAL DURATION: %s   log r2 = %.4f   reg = %.4g x unreg^%.4f"
          % (best["duration"], best["log_r2"], best["log_coeff_a"],
             best["log_exponent_b"]))
    print("   scatter about the fit: x/%.3f (1 sigma, multiplicative)"
          % best["log_se_factor"])

    # ---- cross-check against the flow frequency study's unregulated record --
    if os.path.isfile(MASSBAL_CSV):
        mb = pd.read_csv(MASSBAL_CSV)[["WY", "One_day", "Three_Day", "Five_Day"]]
        mb = mb.rename(columns={"One_day": "mb_One_day",
                                "Three_Day": "mb_Three_Day",
                                "Five_Day": "mb_Five_Day"})
        compare = data.merge(mb, on="WY", how="left")
        if os.path.isfile(WY_PEAKS_CSV):
            wp = pd.read_csv(WY_PEAKS_CSV)[["WY", "unreg_peak_1hr",
                                            "unreg_cov_at_reg_peak"]]
            wp.loc[wp["unreg_cov_at_reg_peak"] < 1.0, "unreg_peak_1hr"] = np.nan
            compare = compare.merge(
                wp.rename(columns={"unreg_peak_1hr": "mb_peak_1hr"})[["WY", "mb_peak_1hr"]],
                on="WY", how="left")
        compare.to_csv(os.path.join(DIAG_DIR, "unreg_source_comparison.csv"),
                       index=False, float_format="%.2f")
        plot_source_comparison(compare, stem)

        print("\nUNREGULATED SOURCE CROSS-CHECK  (ResSim / mass balance)")
        for a_col, b_col, label in [("unreg_Peak_1-hr", "mb_peak_1hr", "Peak (1-hr)"),
                                    ("unreg_1-Day", "mb_One_day", "1-Day"),
                                    ("unreg_3-Day", "mb_Three_Day", "3-Day"),
                                    ("unreg_5-Day", "mb_Five_Day", "5-Day")]:
            if a_col not in compare or b_col not in compare:
                continue
            sub = compare.dropna(subset=[a_col, b_col])
            if len(sub) == 0:
                print("   %-12s no overlapping years" % label)
                continue
            ratio = sub[a_col] / sub[b_col]
            print("   %-12s n=%2d  median ratio %.3f  range %.3f..%.3f"
                  % (label, len(sub), ratio.median(), ratio.min(), ratio.max()))
        print("   Ratios far from 1.0 mean the two unregulated estimates disagree;")
        print("   see %s/unreg_source_comparison.png" % DIAG_DIR)

    print("-" * 78)
    print("Dataset  : %s/critical_duration_adjusted_dataset.csv" % DIAG_DIR)
    print("Fits     : %s/critical_duration_adjusted_fits.csv" % OUT_DIR)
    print("Plots    : %s_scatter.png, %s_summary.png" % (stem, stem))


main()
