"""
Critical_Duration_Correlation.py

Critical duration analysis for Castle Rock: which unregulated flow duration
best explains the observed REGULATED peak?

The idea is that the regulated peak at Castle Rock is controlled by whatever
unregulated inflow duration the Mossyrock/Riffe system cannot absorb. If the
regulated peak tracks the 1-day unregulated flow, the system is peak-limited;
if it tracks the 5-day, the system is volume-limited and the reservoir fills
before the flood ends. The duration with the strongest correlation is the
critical duration.

INPUTS (all read from the CAS_Unreg_FF project -- nothing is duplicated here)

    output/unreg_durations_massbalance.csv
        Unregulated 1/3/5-day water-year maxima from the daily mass balance
        (Castle Rock daily + MOS daily storage change). WY1974-2026.

    data/CastleRock_USGS_peaks.csv
        USGS reported annual peak (the regulated peak). This is the default
        REG_PEAK_SOURCE.

    output/wy_peak_records.csv
        Supplies the unregulated INSTANTANEOUS (1-hr) peak from the hourly
        holdout, included as a zero-day duration for reference, and the
        alternate regulated peak definition (1-hr max of the USGS hourly
        record) if REG_PEAK_SOURCE is switched to "hourly".

SCREENS
    - WY >= FIRST_REG_WY (post-Mossyrock only; no pre-dam years).
    - A water year is only used if the unregulated duration exists. Years
      without unregulated volumes drop out of the correlation by construction.
    - The unregulated 1-hr peak is only used where unreg_cov_at_reg_peak == 1,
      i.e. the holdout actually covered the regulated peak. Low-coverage years
      return a bogus (too low) unregulated peak and would corrupt the fit.
    - EXCLUDE_WYS drops specific water years. WY1980 is excluded by default:
      the 97,000 cfs USGS peak that year is the 18 May 1980 Mount St. Helens
      lahar, not a rainfall-runoff peak, and there is no unregulated volume
      that corresponds to it.

CAVEAT worth keeping in mind when reading the numbers: the unregulated record
is built as regulated + storage change, so the two variables are not
independent by construction. The correlations below are still the right way to
rank durations against each other, but their absolute magnitude should not be
read as skill in a predictive sense.

OUTPUTS (../output and ../diagnostics)
    critical_duration_dataset.csv       joined + screened data, one row per WY
    critical_duration_correlations.csv  fit statistics, one row per duration
    critical_duration_residuals.csv     per-WY residuals for the best fit
    critical_duration_scatter.png       scatter panel with fits, one per duration
    critical_duration_summary.png       r-squared vs duration
"""

import os

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

###############################################################################
# CONFIGURATION

REPO_ROOT = r"C:\Projects\Claude"          # <-- set to your local repo path

SOURCE_PROJECT = os.path.join(REPO_ROOT, "CAS_Unreg_FF")
PROJECT_DIR = os.path.join(REPO_ROOT, "CAS_Reg_Unreg")

durations_csv = os.path.join(SOURCE_PROJECT, "output", "unreg_durations_massbalance.csv")
usgs_peaks_csv = os.path.join(SOURCE_PROJECT, "data", "CastleRock_USGS_peaks.csv")
wy_peaks_csv = os.path.join(SOURCE_PROJECT, "output", "wy_peak_records.csv")

output_dir = os.path.join(PROJECT_DIR, "output")
diag_dir = os.path.join(PROJECT_DIR, "diagnostics")

# "usgs"   -> USGS reported annual peak from CastleRock_USGS_peaks.csv
# "hourly" -> 1-hr max of the USGS hourly record from wy_peak_records.csv
REG_PEAK_SOURCE = "usgs"

FIRST_REG_WY = 1974           # Mossyrock/Riffe regulation begins
EXCLUDE_WYS = [1980]          # see docstring -- Mount St. Helens lahar
INCLUDE_UNREG_INSTANT_PEAK = True   # add the 1-hr unreg peak as a duration
MIN_COV_AT_REG_PEAK = 1.0     # holdout coverage screen for the 1-hr unreg peak

# Column name -> label, in increasing duration order
DURATION_COLS = [
    ("unreg_peak_1hr", "Unreg Peak (1-hr)"),
    ("One_day", "Unreg 1-Day"),
    ("Three_Day", "Unreg 3-Day"),
    ("Five_Day", "Unreg 5-Day"),
]
DURATION_DAYS = {"unreg_peak_1hr": 1.0 / 24.0, "One_day": 1.0, "Three_Day": 3.0, "Five_Day": 5.0}

###############################################################################


def ensure_dirs():
    for d in [output_dir, diag_dir]:
        if not os.path.isdir(d):
            os.makedirs(d)


def build_dataset():
    """Join unregulated durations to the regulated peak and apply screens."""
    dur = pd.read_csv(durations_csv)
    dur = dur[["WY", "One_day", "Three_Day", "Five_Day",
               "days_missing", "flood_season_missing", "Source"]].copy()

    wyp = pd.read_csv(wy_peaks_csv)
    wyp = wyp[["WY", "reg_peak_1hr", "unreg_peak_1hr",
               "unreg_cov_at_reg_peak", "peak_offset_hrs"]].copy()

    usgs = pd.read_csv(usgs_peaks_csv)
    usgs = usgs[["WY", "Peak_cfs", "Peak_Date"]].rename(
        columns={"Peak_cfs": "reg_peak_usgs", "Peak_Date": "reg_peak_usgs_date"})

    df = dur.merge(usgs, on="WY", how="left").merge(wyp, on="WY", how="left")

    # Screen the hourly unregulated peak -- low holdout coverage at the
    # regulated peak means the reported unreg peak is not the real one.
    bad_cov = df["unreg_cov_at_reg_peak"].fillna(-1.0) < MIN_COV_AT_REG_PEAK
    df["unreg_peak_screened_out"] = bad_cov & df["unreg_peak_1hr"].notna()
    df.loc[bad_cov, "unreg_peak_1hr"] = np.nan

    if REG_PEAK_SOURCE == "usgs":
        df["reg_peak"] = df["reg_peak_usgs"]
        df["reg_peak_source"] = np.where(df["reg_peak_usgs"].notna(), "usgs_peak_record", "")
    else:
        df["reg_peak"] = df["reg_peak_1hr"]
        df["reg_peak_source"] = np.where(df["reg_peak_1hr"].notna(), "hourly_1hr_max", "")

    df["used"] = True
    df.loc[df["WY"] < FIRST_REG_WY, "used"] = False
    df.loc[df["WY"].isin(EXCLUDE_WYS), "used"] = False
    df.loc[df["reg_peak"].isna(), "used"] = False
    df.loc[df["Three_Day"].isna(), "used"] = False

    df["exclude_reason"] = ""
    df.loc[df["WY"] < FIRST_REG_WY, "exclude_reason"] = "pre-regulation"
    df.loc[df["WY"].isin(EXCLUDE_WYS), "exclude_reason"] = "manual EXCLUDE_WYS"
    df.loc[df["reg_peak"].isna(), "exclude_reason"] = "no regulated peak"
    df.loc[df["Three_Day"].isna(), "exclude_reason"] = "no unregulated volume"

    return df.sort_values("WY").reset_index(drop=True)


def fit_one(x, y):
    """OLS plus correlation measures for one predictor. Returns a dict."""
    n = len(x)
    slope, intercept, r, p, stderr = stats.linregress(x, y)
    resid = y - (slope * x + intercept)
    se_est = float(np.sqrt(np.sum(resid ** 2) / (n - 2))) if n > 2 else np.nan
    rho, rho_p = stats.spearmanr(x, y)
    tau, tau_p = stats.kendalltau(x, y)
    return {
        "n": n,
        "pearson_r": r,
        "r2": r ** 2,
        "p_value": p,
        "spearman_rho": rho,
        "kendall_tau": tau,
        "slope": slope,
        "intercept": intercept,
        "slope_stderr": stderr,
        "se_estimate_cfs": se_est,
    }


def fit_log(x, y):
    """Power-law fit in log10 space: reg = a * dur ** b."""
    lx = np.log10(x)
    ly = np.log10(y)
    slope, intercept, r, p, stderr = stats.linregress(lx, ly)
    resid = ly - (slope * lx + intercept)
    n = len(lx)
    se_log = float(np.sqrt(np.sum(resid ** 2) / (n - 2))) if n > 2 else np.nan
    return {
        "log_r2": r ** 2,
        "log_exponent_b": slope,
        "log_coeff_a": 10.0 ** intercept,
        "log_se": se_log,
    }


def active_columns():
    cols = []
    for col, label in DURATION_COLS:
        if col == "unreg_peak_1hr" and not INCLUDE_UNREG_INSTANT_PEAK:
            continue
        cols.append((col, label))
    return cols


def run_correlations(df, subset):
    """One row of statistics per candidate duration.

    subset = "all_available": each duration uses every water year it has.
    subset = "common_years":  every duration is fit on the identical set of
                              water years, so the r-squared values are
                              directly comparable. The 1-hr unregulated peak
                              has far fewer years than the daily durations, so
                              the all_available ranking mixes sample sizes.
    """
    use = df[df["used"]].copy()
    cols = active_columns()
    if subset == "common_years":
        needed = ["reg_peak"] + [c for c, _ in cols]
        use = use.dropna(subset=needed)

    rows = []
    for col, label in cols:
        sub = use[["WY", col, "reg_peak"]].dropna()
        if len(sub) < 3:
            continue
        x = sub[col].to_numpy(dtype=float)
        y = sub["reg_peak"].to_numpy(dtype=float)
        row = {"subset": subset, "duration": label, "column": col,
               "duration_days": DURATION_DAYS[col],
               "first_wy": int(sub["WY"].min()), "last_wy": int(sub["WY"].max())}
        row.update(fit_one(x, y))
        row.update(fit_log(x, y))
        row["mean_reg_over_unreg"] = float(np.mean(y / x))
        rows.append(row)
    stats_df = pd.DataFrame(rows)
    if stats_df.empty:
        return stats_df
    return stats_df.sort_values("duration_days").reset_index(drop=True)


def residual_table(df, col):
    """Per-WY residuals and standardized residuals for one duration."""
    sub = df[df["used"]][["WY", col, "reg_peak"]].dropna().copy()
    x = sub[col].to_numpy(dtype=float)
    y = sub["reg_peak"].to_numpy(dtype=float)
    slope, intercept, _, _, _ = stats.linregress(x, y)
    sub["predicted_reg_peak"] = slope * x + intercept
    sub["residual_cfs"] = sub["reg_peak"] - sub["predicted_reg_peak"]
    sd = sub["residual_cfs"].std(ddof=2)
    sub["standardized_residual"] = sub["residual_cfs"] / sd if sd > 0 else np.nan
    sub["outlier_flag"] = np.abs(sub["standardized_residual"]) > 2.0
    return sub.sort_values("standardized_residual").reset_index(drop=True)


def plot_scatter(df, stats_df, best, path):
    use = df[df["used"]]
    cols = [(c, l) for c, l in DURATION_COLS
            if c in set(stats_df["column"])]
    ncol = len(cols)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 4.4), sharey=True)
    if ncol == 1:
        axes = [axes]

    for ax, (col, label) in zip(axes, cols):
        sub = use[["WY", col, "reg_peak"]].dropna()
        srow = stats_df[stats_df["column"] == col].iloc[0]
        x = sub[col].to_numpy(dtype=float)
        y = sub["reg_peak"].to_numpy(dtype=float)
        ax.scatter(x, y, s=26, color="#1f4e79", alpha=0.8, zorder=3)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, srow["slope"] * xs + srow["intercept"],
                color="#c00000", lw=1.6, zorder=2)
        ax.plot(xs, xs, color="0.6", lw=1.0, ls="--", zorder=1)
        title = label if col != best else label + "  (best)"
        ax.set_title(title, fontsize=11,
                     fontweight="bold" if col == best else "normal")
        ax.set_xlabel("Unregulated flow (cfs)")
        ax.grid(True, ls=":", lw=0.6, color="0.8")
        txt = ("n = %d\nr$^2$ = %.3f\n$\\rho$ = %.3f\nSE = %s cfs"
               % (srow["n"], srow["r2"], srow["spearman_rho"],
                  format(int(round(srow["se_estimate_cfs"])), ",")))
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
                fontsize=8.5, bbox=dict(fc="white", ec="0.7", alpha=0.9, pad=4))
    axes[0].set_ylabel("Regulated peak at Castle Rock (cfs)")
    fig.suptitle("Critical duration: regulated peak vs unregulated duration "
                 "(WY%d-%d, dashed line is 1:1)"
                 % (int(use["WY"].min()), int(use["WY"].max())), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_summary(stats_df, path):
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    labels = stats_df["duration"].tolist()
    idx = np.arange(len(labels))
    w = 0.27
    ax.bar(idx - w, stats_df["r2"], w, label="Pearson r$^2$ (linear)", color="#1f4e79")
    ax.bar(idx, stats_df["log_r2"], w, label="Pearson r$^2$ (log-log)", color="#4472c4")
    ax.bar(idx + w, stats_df["spearman_rho"] ** 2, w, label="Spearman $\\rho^2$", color="#a6a6a6")
    ax.set_xticks(idx)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Explained variance")
    ax.set_ylim(0, 1.0)
    ax.grid(True, axis="y", ls=":", lw=0.6, color="0.8")
    ax.set_axisbelow(True)
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("Correlation of regulated peak with unregulated duration", fontsize=12)
    for i, v in enumerate(stats_df["r2"]):
        ax.text(i - w, v + 0.015, "%.3f" % v, ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def print_stats_table(stats_df):
    show = ["duration", "n", "pearson_r", "r2", "log_r2", "spearman_rho",
            "kendall_tau", "slope", "intercept", "se_estimate_cfs"]
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(stats_df[show].to_string(
            index=False,
            formatters={"pearson_r": "{:.4f}".format, "r2": "{:.4f}".format,
                        "log_r2": "{:.4f}".format, "spearman_rho": "{:.4f}".format,
                        "kendall_tau": "{:.4f}".format, "slope": "{:.4f}".format,
                        "intercept": "{:,.0f}".format,
                        "se_estimate_cfs": "{:,.0f}".format}))


def main():
    ensure_dirs()

    df = build_dataset()
    stats_all = run_correlations(df, "all_available")
    stats_common = run_correlations(df, "common_years")
    if stats_all.empty:
        print("No duration had enough paired years to fit. Check the screens.")
        return
    stats_df = pd.concat([stats_all, stats_common], ignore_index=True)

    # Rank only the daily durations -- those are the critical-duration
    # candidates and they share the same water years. The 1-hr unregulated
    # peak is carried for reference but sits on a much shorter record.
    daily = stats_all[stats_all["column"] != "unreg_peak_1hr"]
    ranked = daily if not daily.empty else stats_all
    best_col = ranked.loc[ranked["r2"].idxmax(), "column"]
    best_label = ranked.loc[ranked["r2"].idxmax(), "duration"]
    resid = residual_table(df, best_col)

    ds_path = os.path.join(output_dir, "critical_duration_dataset.csv")
    st_path = os.path.join(output_dir, "critical_duration_correlations.csv")
    rs_path = os.path.join(diag_dir, "critical_duration_residuals.csv")
    sc_path = os.path.join(output_dir, "critical_duration_scatter.png")
    sm_path = os.path.join(output_dir, "critical_duration_summary.png")

    df.to_csv(ds_path, index=False)
    stats_df.to_csv(st_path, index=False)
    resid.to_csv(rs_path, index=False)
    plot_scatter(df, stats_all, best_col, sc_path)
    plot_summary(stats_all, sm_path)

    used = df[df["used"]]
    dropped = df[~df["used"]]

    print("=" * 78)
    print("CRITICAL DURATION ANALYSIS -- Castle Rock (USGS 14243000)")
    print("=" * 78)
    print("Regulated peak source : %s" % REG_PEAK_SOURCE)
    print("Water years used      : %d  (WY%d-%d)"
          % (len(used), used["WY"].min(), used["WY"].max()))
    if len(dropped):
        print("Water years dropped   : %d" % len(dropped))
        for _, r in dropped.iterrows():
            print("    WY%-5d %s" % (r["WY"], r["exclude_reason"]))
    n_screen = int(df["unreg_peak_screened_out"].sum())
    if n_screen:
        print("Unreg 1-hr peaks screened out for low holdout coverage: %d" % n_screen)
    print()

    print("ALL AVAILABLE YEARS (each duration uses every year it has)")
    print_stats_table(stats_all)
    if not stats_common.empty:
        n_common = int(stats_common["n"].iloc[0])
        print()
        print("COMMON YEARS ONLY (n = %d, identical sample for every duration)" % n_common)
        print_stats_table(stats_common)
    print()
    print("Best-correlated daily duration (linear r^2): %s" % best_label)
    r2_span = ranked["r2"].max() - ranked["r2"].min()
    print("Spread in r^2 across the daily durations: %.4f" % r2_span)
    if r2_span < 0.05:
        print("  -> The durations are nearly indistinguishable. Treat the ranking")
        print("     as weak evidence and lean on the physical argument.")
    print()

    out = resid[resid["outlier_flag"]]
    if len(out):
        print("Water years with |standardized residual| > 2 on the %s fit:" % best_label)
        for _, r in out.iterrows():
            print("    WY%-5d resid %+10s cfs   z = %+.2f"
                  % (r["WY"], format(int(round(r["residual_cfs"])), ","),
                     r["standardized_residual"]))
    else:
        print("No |standardized residual| > 2 on the %s fit." % best_label)
    print()
    print("Wrote:")
    for p in [ds_path, st_path, rs_path, sc_path, sm_path]:
        print("    %s" % p)


main()
