"""
PeakDiff_Storage_Regression.py

Regress the water-year peak difference (REG minus UNREG, 1-hr peaks from
WY_Peak_Records.py) against the daily storage change at Mossyrock, and
apply the winning regression to estimate an unregulated peak for WYs
that have a good regulated peak but no usable hourly holdout.

Storage-change predictors are computed from the daily MOS elevation
record (//MOS/ELEV//1DAY/USGS/ -- much cleaner than the hourly
telemetry), converted to storage with the official 2014 rating. For each WY, within a window around the
regulated peak date, several candidate metrics are computed:

    dS_1day   max 1-day storage increase
    dS_2day   max 2-day storage increase
    dS_3day   max 3-day storage increase
    dS_4day   max 4-day storage increase  (== greatest change in any
                                           4-day period in the window)

Each is reported in ac-ft and as an equivalent mean flow in cfs over
its window length. A separate OLS fit is run for each candidate
predictor; slope/intercept/R^2/SE are tabulated so the best-performing
window can be picked. WYs whose reg/unreg peaks landed far apart in
time (see MAX_PEAK_OFFSET_HRS) are excluded from the fit -- those pairs
are probably different storms.

Application: for WYs with a regulated peak (from wy_peak_records.csv,
or the USGS peak-flow CSV as a fallback) but no unreg peak, the chosen
regression predicts REG-minus-UNREG from the WY's dS metric, and
    unreg_peak_est = reg_peak - predicted_diff
(predicted_diff is normally negative during floods -- storage capture
makes unreg exceed reg -- so the estimate lands above the reg peak).

Outputs (../output and ../output/diagnostics):
    peakdiff_storage_regressions.csv   fit stats for every candidate
    peakdiff_regression_pairs.csv      the fit dataset (audit trail)
    unreg_peak_estimates.csv           corrected peaks for gap WYs
    peakdiff_storage_regression.png    scatter panel, one per candidate
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

###############################################################################
# CONFIGURATION

# Repo root derived from this file: <repo>/CAS_Unreg_FF/src/<script>.py
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
root_dir = os.path.join(PROJECT_DIR, "data")
peaks_dir = os.path.join(PROJECT_DIR, "output")
output_dir = os.path.join(PROJECT_DIR, "output")
diag_dir = os.path.join(PROJECT_DIR, "output", "diagnostics")

sys.path.insert(0, REPO_ROOT)
MODULES_DIR = os.path.join(REPO_ROOT, "Modules")
sys.path.insert(0, MODULES_DIR)
from utilsDSS import HecDss  # noqa: E402
# Official 2014 elev<->stor rating + interp/extrap live in the holdout script
from Build_Hourly_Holdout_Unreg import elev_to_stor  # noqa: E402

DSS_OBS = os.path.join(root_dir, "obsData.dss")
PEAKS_CSV = os.path.join(peaks_dir, "wy_peak_records.csv")

# Daily MOS elevation for the storage-change predictors: the daily
# USGS record -- a separate, cleaner record, NOT derived from the
# hourly CWMS-CLEAN telemetry (much of the post-1974 hourly record is
# not usable at hourly resolution; CWMS-CLEAN is used only in the
# holdout workflow, gated day-by-day by the STOR-COUNT screen and
# manual overrides). The hourly-derived option remains available via
# the flag for comparison runs only.
DAILY_FROM_HOURLY_CLEAN = False
PATH_MOS_ELEV_HOURLY_CLEAN = "//MOS/ELEV//1HOUR/CWMS-CLEAN/"
PATH_MOS_ELEV_DAILY = "//MOS/ELEV//1DAY/USGS/"

# Optional fallback source of regulated peaks for gap WYs where the USGS
# hourly record is not good either (USGS instantaneous annual peaks).
# Set to None to disable.
USGS_PEAKS_CSV = os.path.join(PROJECT_DIR, "data",
                              "CastleRock_USGS_peaks.csv")

OUT_FITS_CSV = os.path.join(output_dir, "peakdiff_storage_regressions.csv")
OUT_PAIRS_CSV = os.path.join(output_dir, "peakdiff_regression_pairs.csv")
OUT_EST_CSV = os.path.join(output_dir, "unreg_peak_estimates.csv")
OUT_PNG = os.path.join(diag_dir, "peakdiff_storage_regression.png")

WINDOW_DAYS = 7            # +/- days around the reg peak date to search for dS
# Gap test for WYs that DO have an hourly unreg peak: the peak is
# unusable when the unregulated record has no data covering the basin's
# biggest event, as dated by the regulated annual peak
# (unreg_cov_at_reg_peak from WY_Peak_Records). Same threshold and same
# column Write_SSP_Record uses, so the two scripts always agree.
MIN_EVENT_COVERAGE = 0.9
REG_FILL_FIRST_WY = 1969   # regression fill applies to the regulated era only
DS_WINDOWS = [1, 2, 3, 4]  # storage-change window lengths, days
MAX_PEAK_OFFSET_HRS = 72   # exclude WYs whose reg/unreg peaks are farther
                           # apart than this (likely different storms)
CFS_DAYS_PER_AF = 43560.0 / 86400.0  # 1 ac-ft over 1 day = 0.504 cfs

# Which fitted predictor to use for the applied correction. One of
# "dS_1day_cfs" ... "dS_4day_cfs", or None to auto-pick the highest R^2.
# ADOPTED: dS_2day (best-performing window, Jul 2026).
APPLY_PREDICTOR = "dS_2day_cfs"

###############################################################################
# FUNCTION DEFINITIONS


def read_daily_stor():
    """Read the daily MOS elevation (either the daily mean of the
    hand-cleaned hourly CWMS-CLEAN record, or a designated daily record)
    and convert to storage (ac-ft)."""
    path = PATH_MOS_ELEV_HOURLY_CLEAN if DAILY_FROM_HOURLY_CLEAN \
        else PATH_MOS_ELEV_DAILY
    if not path:
        raise RuntimeError(
            "PATH_MOS_ELEV_DAILY is empty -- set it, or set "
            "DAILY_FROM_HOURLY_CLEAN = True to derive daily means from "
            "the CWMS-CLEAN hourly record.")
    dss = HecDss.open(DSS_OBS)
    try:
        df = dss.readDF(path)
    finally:
        dss.close()
    if df.empty:
        raise RuntimeError(f"No data for {path} in {DSS_OBS}")
    elev = df["value"]
    elev = elev.mask(elev <= -900.0)
    elev = elev[~elev.index.duplicated(keep="last")].sort_index()
    elev = elev.resample("1D").mean()  # hourly -> daily mean; daily -> grid
    stor = elev_to_stor(elev)
    stor[elev.isna()] = np.nan
    print(f"  daily STOR: {stor.notna().sum()} valid days, "
          f"{stor.index[0].date()} -> {stor.index[-1].date()}")
    return stor


def ds_metrics(stor, peak_date):
    """Storage-change metrics near one event. For each N in DS_WINDOWS,
    the max of stor[t] - stor[t-N] over the +/-WINDOW_DAYS window around
    peak_date. Returns dict of {name: value} in ac-ft and mean-cfs."""
    if pd.isna(peak_date):
        return {}
    d0 = pd.Timestamp(peak_date).normalize()
    out = {}
    for n in DS_WINDOWS:
        win = stor.loc[d0 - pd.Timedelta(days=WINDOW_DAYS):
                       d0 + pd.Timedelta(days=WINDOW_DAYS)]
        diff = win - win.shift(n)
        if diff.notna().sum() == 0:
            out[f"dS_{n}day_af"] = np.nan
            out[f"dS_{n}day_cfs"] = np.nan
            continue
        v = diff.max()
        out[f"dS_{n}day_af"] = v
        out[f"dS_{n}day_cfs"] = v * CFS_DAYS_PER_AF / n  # mean cfs over N days
    return out


def fit_ols(x, y):
    """Simple OLS y = a*x + b. Returns dict with slope, intercept, r2,
    se (std of residuals), n."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 3:
        return {"slope": np.nan, "intercept": np.nan, "r2": np.nan,
                "se": np.nan, "n": n}
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    resid = y - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    se = float(np.sqrt(ss_res / (n - 2)))
    return {"slope": a, "intercept": b, "r2": r2, "se": se, "n": n}


def load_usgs_peaks():
    """USGS instantaneous annual peaks: DataFrame indexed by WY with
    Peak_cfs and Peak_Date. The Peak_Date matters -- it centers the dS
    search window for WYs that never enter the hourly peak table."""
    if USGS_PEAKS_CSV is None or not os.path.exists(USGS_PEAKS_CSV):
        return pd.DataFrame(
            columns=["Peak_cfs", "Peak_Date"]).rename_axis("WY")
    df = pd.read_csv(USGS_PEAKS_CSV, parse_dates=["Peak_Date"])
    return df.set_index("WY")[["Peak_cfs", "Peak_Date"]]


def plot_panels(pairs, fits):
    """One scatter per candidate predictor with the fitted line."""
    preds = [f"dS_{n}day_cfs" for n in DS_WINDOWS]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharey=True)
    for ax, p in zip(axes.ravel(), preds):
        x = pairs[p].values
        y = pairs["reg_minus_unreg_1hr"].values
        ax.scatter(x, y, s=25, edgecolor="k", linewidth=0.4, zorder=3)
        for wy, xi, yi in zip(pairs.index, x, y):
            if np.isfinite(xi) and np.isfinite(yi):
                ax.annotate(str(wy), (xi, yi), fontsize=6.5,
                            xytext=(3, 3), textcoords="offset points")
        f = fits.loc[p]
        if np.isfinite(f["slope"]):
            xr = np.linspace(np.nanmin(x), np.nanmax(x), 10)
            ax.plot(xr, f["slope"] * xr + f["intercept"], "r-", lw=1.2)
        ax.set_title(f"{p}: R\u00b2={f['r2']:.3f}, SE={f['se']:,.0f} cfs, "
                     f"n={int(f['n'])}", fontsize=9)
        ax.set_xlabel("max storage change (mean cfs over window)")
        ax.set_ylabel("REG - UNREG peak (cfs)")
        ax.grid(alpha=0.3)
    fig.suptitle("Castle Rock peak difference vs MOS daily storage change",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"  wrote {OUT_PNG}")


###############################################################################
# MAIN


def main():
    print(f"Reading WY peak table {PEAKS_CSV}")
    peaks = pd.read_csv(PEAKS_CSV, index_col="WY",
                        parse_dates=["reg_peak_1hr_time",
                                     "unreg_peak_1hr_time"])
    print(f"Reading clean daily MOS elevation from {DSS_OBS}")
    stor = read_daily_stor()

    usgs = load_usgs_peaks()

    # --- storage-change metrics for every WY with ANY peak date:
    # hourly reg peak time where the hourly table has one, else the
    # USGS peak-record date (regulated era only) ---
    all_wys = sorted(set(peaks.index)
                     | set(usgs.index[usgs.index >= REG_FILL_FIRST_WY]))
    date_map = {}
    for wy in all_wys:
        t = peaks["reg_peak_1hr_time"].get(wy, pd.NaT) \
            if wy in peaks.index else pd.NaT
        if pd.isna(t) and wy in usgs.index:
            t = usgs.loc[wy, "Peak_Date"]
        date_map[wy] = t
    metric_rows = {wy: ds_metrics(stor, date_map[wy]) for wy in all_wys}
    metrics = pd.DataFrame.from_dict(metric_rows, orient="index")
    metrics.index.name = "WY"
    table = peaks.join(metrics)  # fit uses hourly-table years only

    # --- fit dataset: WYs with both peaks, same-storm screen ---
    fitset = table.dropna(subset=["reg_minus_unreg_1hr"])
    off = fitset["peak_offset_hrs"].abs()
    dropped = fitset[off > MAX_PEAK_OFFSET_HRS]
    if len(dropped):
        print(f"  excluding {len(dropped)} WY(s) with peak offset > "
              f"{MAX_PEAK_OFFSET_HRS} hrs: {list(dropped.index)}")
    fitset = fitset[off <= MAX_PEAK_OFFSET_HRS]
    fitset.to_csv(OUT_PAIRS_CSV)
    print(f"  fit dataset: {len(fitset)} WYs -> {OUT_PAIRS_CSV}")

    # --- one OLS per candidate predictor ---
    fits = {}
    for n in DS_WINDOWS:
        p = f"dS_{n}day_cfs"
        fits[p] = fit_ols(fitset[p].values,
                          fitset["reg_minus_unreg_1hr"].values)
    fits = pd.DataFrame.from_dict(fits, orient="index")
    fits.index.name = "predictor"
    fits.to_csv(OUT_FITS_CSV)
    print("\nRegression candidates:")
    print(fits.to_string(float_format=lambda v: f"{v:,.3f}"))

    plot_panels(fitset, fits)

    # --- apply the chosen regression to gap WYs ---
    pred_name = APPLY_PREDICTOR or fits["r2"].idxmax()
    f = fits.loc[pred_name]
    print(f"\nApplying predictor {pred_name} "
          f"(R\u00b2={f['r2']:.3f}) to gap WYs")

    # gap universe (regulated era, WY >= REG_FILL_FIRST_WY):
    #   no_unreg_peak hourly table row exists but no unreg peak
    #   missed_event  unreg peak exists but the unreg record has no data
    #                 covering the regulated annual peak event -- the
    #                 flood fell in a gap in the hourly record
    #   no_hourly     WY only in the USGS peak record (hourly record
    #                 hasn't started / has no row) -- e.g. WY1974-1991
    evcov = table["unreg_cov_at_reg_peak"] \
        if "unreg_cov_at_reg_peak" in table.columns \
        else pd.Series(np.nan, index=table.index)
    no_unreg = set(table.index[table["unreg_peak_1hr"].isna()])
    missed = set(table.index[table["unreg_peak_1hr"].notna()
                             & ~(evcov >= MIN_EVENT_COVERAGE)])
    csv_only = set(usgs.index[usgs.index >= REG_FILL_FIRST_WY]) \
        - set(table.index)
    gap_wys = sorted(y for y in (no_unreg | missed | csv_only)
                     if y >= REG_FILL_FIRST_WY)
    print(f"  gap WYs: {len(no_unreg)} no-unreg-peak, {len(missed)} "
          f"missed-event(event coverage < {MIN_EVENT_COVERAGE}), "
          f"{len(csv_only)} USGS-peak-record-only")

    est_rows = []
    for wy in gap_wys:
        if wy in table.index \
                and np.isfinite(table["reg_peak_1hr"].get(wy, np.nan)):
            reg_peak = float(table.loc[wy, "reg_peak_1hr"])
            source = "hourly_1hr_max"
        elif wy in usgs.index \
                and np.isfinite(usgs.loc[wy, "Peak_cfs"]):
            reg_peak = float(usgs.loc[wy, "Peak_cfs"])
            source = "usgs_peak_record"
        else:
            continue
        x = metrics[pred_name].get(wy, np.nan)
        if not (np.isfinite(reg_peak) and np.isfinite(x)):
            print(f"    WY{wy}: reg peak available but no dS predictor "
                  "(daily elevation window empty) -- not filled")
            continue
        diff_pred = f["slope"] * x + f["intercept"]
        reason = ("no_hourly" if wy in csv_only else
                  "missed_event" if wy in missed else "no_unreg_peak")
        est_rows.append({
            "WY": wy,
            "reg_peak": reg_peak,
            "reg_peak_source": source,
            "gap_reason": reason,
            pred_name: x,
            "predicted_reg_minus_unreg": diff_pred,
            "unreg_peak_est": reg_peak - diff_pred,
            "se_cfs": f["se"],
        })
    est = pd.DataFrame(est_rows)
    est.to_csv(OUT_EST_CSV, index=False)
    print(f"  {len(est)} gap-WY estimates -> {OUT_EST_CSV}")


if __name__ == "__main__":
    main()
