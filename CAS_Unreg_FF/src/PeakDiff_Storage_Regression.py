"""
PeakDiff_Storage_Regression.py

Regress the water-year peak difference (REG minus UNREG, 1-hr peaks from
WY_Peak_Records.py) against the daily storage change at Mossyrock, and
apply the winning regression to estimate an unregulated peak for WYs
that have a good regulated peak but no usable hourly holdout.

Storage-change predictors are computed from a CLEAN DAILY MOS elevation
record (much cleaner than the hourly record), converted to storage with
the official 2014 rating. For each WY, within a window around the
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

Outputs (../output and ../diagnostics):
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

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path
USE_REFERENCE_DATA = True                  # True = ref_data sample run

PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Unreg_FF")

if USE_REFERENCE_DATA:
    root_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_in")
    peaks_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_out")
    output_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_out")
    diag_dir = os.path.join(PROJECT_DIR, "ref_data", "ref_out")
else:
    root_dir = os.path.join(PROJECT_DIR, "data")
    peaks_dir = os.path.join(PROJECT_DIR, "output")
    output_dir = os.path.join(PROJECT_DIR, "output")
    diag_dir = os.path.join(PROJECT_DIR, "diagnostics")

sys.path.insert(0, REPO_ROOT)
UTILS_DIR = os.path.join(PROJECT_DIR, "src", "Cowlitz_Unreg", "Cowlitz")
sys.path.insert(0, UTILS_DIR)
from utilsDSS import HecDss  # noqa: E402
# Official 2014 elev<->stor rating + interp/extrap live in the holdout script
from Build_Hourly_Holdout_Unreg import elev_to_stor  # noqa: E402

DSS_OBS = os.path.join(root_dir, "obsData.dss")
PEAKS_CSV = os.path.join(peaks_dir, "wy_peak_records.csv")

# Daily MOS elevation source for the storage-change predictors.
#   Default: daily mean of the hand-cleaned hourly record (CWMS-CLEAN) --
#   the record actually maintained under the current methodology.
#   If a separate, cleaner pre-computed DAILY record is designated later,
#   set PATH_MOS_ELEV_DAILY and flip DAILY_FROM_HOURLY_CLEAN to False.
#   (The legacy //MOS/ELEV-FOREBAY//1DAY/IRVZZAZD_CLEANED/ record is
#   retired and slated for deletion -- do not point here.)
DAILY_FROM_HOURLY_CLEAN = True
PATH_MOS_ELEV_HOURLY_CLEAN = "//MOS/ELEV//1HOUR/CWMS-CLEAN/"
PATH_MOS_ELEV_DAILY = ""  # only used when DAILY_FROM_HOURLY_CLEAN = False

# Optional fallback source of regulated peaks for gap WYs where the USGS
# hourly record is not good either (USGS instantaneous annual peaks).
# Set to None to disable.
USGS_PEAKS_CSV = os.path.join(REPO_ROOT, "Cowlitz_FF_DataPrep", "data",
                              "CastleRock_USGS_peaks.csv")

OUT_FITS_CSV = os.path.join(output_dir, "peakdiff_storage_regressions.csv")
OUT_PAIRS_CSV = os.path.join(output_dir, "peakdiff_regression_pairs.csv")
OUT_EST_CSV = os.path.join(output_dir, "unreg_peak_estimates.csv")
OUT_PNG = os.path.join(diag_dir, "peakdiff_storage_regression.png")

WINDOW_DAYS = 7            # +/- days around the reg peak date to search for dS
DS_WINDOWS = [1, 2, 3, 4]  # storage-change window lengths, days
MAX_PEAK_OFFSET_HRS = 72   # exclude WYs whose reg/unreg peaks are farther
                           # apart than this (likely different storms)
CFS_DAYS_PER_AF = 43560.0 / 86400.0  # 1 ac-ft over 1 day = 0.504 cfs

# Which fitted predictor to use for the applied correction. One of
# "dS_1day_cfs" ... "dS_4day_cfs", or None to auto-pick the highest R^2.
APPLY_PREDICTOR = None

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
    """Optional USGS instantaneous annual peaks (WY, peak_cfs, date)."""
    if USGS_PEAKS_CSV is None or not os.path.exists(USGS_PEAKS_CSV):
        return pd.DataFrame()
    df = pd.read_csv(USGS_PEAKS_CSV)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


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

    # --- storage-change metrics for every WY with a reg peak date ---
    metric_rows = {}
    for wy, row in peaks.iterrows():
        metric_rows[wy] = ds_metrics(stor, row.get("reg_peak_1hr_time"))
    metrics = pd.DataFrame.from_dict(metric_rows, orient="index")
    metrics.index.name = "WY"
    table = peaks.join(metrics)

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

    gaps = table[table["unreg_peak_1hr"].isna()].copy()
    usgs = load_usgs_peaks()
    est_rows = []
    for wy, row in gaps.iterrows():
        reg_peak = row.get("reg_peak_1hr", np.nan)
        source = "hourly_1hr_max"
        if not np.isfinite(reg_peak) and len(usgs):
            hit = usgs[usgs.iloc[:, 0] == wy]  # first col assumed WY
            if len(hit):
                num_cols = hit.select_dtypes("number").columns
                if len(num_cols) > 1:
                    reg_peak = float(hit.iloc[0][num_cols[1]])
                    source = "usgs_peak_record"
        x = row.get(pred_name, np.nan)
        if not (np.isfinite(reg_peak) and np.isfinite(x)):
            continue
        diff_pred = f["slope"] * x + f["intercept"]
        est_rows.append({
            "WY": wy,
            "reg_peak": reg_peak,
            "reg_peak_source": source,
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
