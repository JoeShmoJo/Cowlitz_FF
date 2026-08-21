#Coweeman_FlowFrequency.py
# -*- coding: utf-8 -*-
"""
Coweeman River unregulated peak flow-frequency curve.

WHY THIS EXISTS
    The Castle Rock coincident-frequency work (see
    CAS_Reg_Unreg/docs/CDID3_Coincident_Frequency_Notes.md) needs a Coweeman
    marginal frequency curve for two of the combination methods --
    #Coincident_PerfectCorrelation.py (same-AEP sum, the CDID3 Phase 1
    precedent) and #Coincident_CorrConditioned.py (correlation-conditioned
    combination). This script builds that curve so those two have something
    real to run against.

WHAT THIS IS, AND WHAT IT IS NOT
    The CDID3 Phase 2a report (2016) built its own Coweeman curve in
    HEC-SSP: Bulletin 17C draft methods, multiple Grubbs-Beck low-outlier
    screening, the Expected Moments Algorithm (EMA) to weight the 1933 and
    1996 historic peaks and the 1985-2006 perception threshold against the
    systematic record, and a regional skew from Cooper (2005). That is real
    statistical machinery -- this repo's OWN unregulated Cowlitz curve goes
    through the same kind of process. See
    CAS_Unreg_FF/src/Frequency_Curves_And_Table.py, which only PARSES an
    HEC-SSP report; the actual fit happens in SSP, not in Python.

    HEC-SSP is not available in this environment. Rather than hand-roll an
    EMA / multiple-Grubbs-Beck implementation from memory -- exactly the
    kind of thing that produces a plausible-looking curve with a subtle,
    hard-to-catch error in it -- this script does something more modest and
    easier to check by hand:
      - Log-Pearson III fit by ordinary method of moments, on the
        SYSTEMATIC record alone (USGS 1950-1984 annual peaks, plus
        Ecology-derived water-year maxima 2007-2020).
      - Station skew blended with the SAME regional skew CDID3 used
        (0.2, MSE 0.112, Cooper 2005), via the standard Bulletin 17B
        weighting formula (a closed-form algebraic formula, not an
        iterative EMA solve).
      - The 1933 and 1996 historic peaks are PLOTTED at a historic-adjusted
        plotting position for context, but are NOT folded into the moment
        estimates -- no historical weighting, no EMA.
      - Confidence limits by bootstrap (resample the systematic record with
        replacement, refit, repeat many times, take percentiles) rather
        than the WRC/Bulletin 17B noncentral-t confidence-limit tables.

    This is an interim approximation, not a substitute for an actual
    HEC-SSP run. It should be expected to disagree with a real EMA fit,
    especially at the low-AEP tail, where the perception threshold and the
    two historic peaks would carry real weight in a full Bulletin 17C
    treatment and carry none here. It exists so the coincident-frequency
    scripts have real numbers to develop and test against now. BEFORE this
    curve goes in the memo: run the same record through HEC-SSP the way the
    Cowlitz curve was, and turn this script into a report parser the way
    Frequency_Curves_And_Table.py is for Cowlitz, rather than a fitter.

LOCATION -- ONE OPEN ITEM
    This curve is left at the gage (USGS 14245000, "Coweeman River near
    Kelso"), with no drainage-area adjustment applied. CDID3 scaled their
    own gage-based curve up by a factor of 1.07 (127 vs 119 sq mi) to reach
    their levee, which sits a short distance further downstream than the
    gage. Whether an equivalent correction is needed here depends on where
    exactly the coincident combination is evaluated relative to the gage,
    which isn't settled yet -- flagged for review rather than guessed at.

DATA
    Systematic annual peaks:
      1950-1984  USGS gage 14245000, from
                 data/coweeman/usgs_peaks_14245000.rdb
      2007-2020  Water-year max of the Ecology 15-minute record at 26C075,
                 from data/coweeman/ecology_26C075_*_FM.txt (the same files
                 #Coweeman_Timing.py fetches and caches)
    The gage was inactive 1985-2006 except for one measured event:
      1996-02-08  11,700 cfs (USGS 14245000, flagged peak_cd=7 in the rdb)
    Historic peak, pre-systematic record, per the CDID3 report citing FEMA
    (2015):
      1933        12,000 cfs
    Neither historic value is independently re-derived here -- both are
    literal numbers transcribed from the CDID3 Phase 2a report text (see
    CAS_Reg_Unreg/docs/CDID3_Coincident_Frequency_Notes.md).

REGIONAL SKEW
    0.2, MSE 0.112 -- USGS regional skew study at the Coweeman gage
    (Cooper 2005), the same value CDID3 cites and uses.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
PEAKS_RDB = r"../data/coweeman/usgs_peaks_14245000.rdb"
ECOLOGY_GLOB = r"../data/coweeman/ecology_26C075_*_FM.txt"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "coweeman_frequency_table.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coweeman_frequency_curve.png")

# Historic peaks -- transcribed from the CDID3 Phase 2a report (2016), not
# independently sourced. Shown on the plot for context, not fit.
HISTORIC_PEAKS = {
    1933: 12000.0,    # FEMA (2015), via the CDID3 report
    1996: 11700.0,    # USGS 14245000, peak_cd=7 (also present in the rdb)
}
HISTORIC_PERIOD_START = 1933     # for the historic-adjusted plotting position

REGIONAL_SKEW = 0.2
REGIONAL_SKEW_MSE = 0.112        # Cooper (2005)

# Matches the AEP grid already used in CAS_Unreg_FF/output/
# CAS_Unreg_frequency_table.csv and CAS_Reg_Unreg/output/
# regulated_frequency_inferred.csv, so the three curves interpolate onto a
# shared grid with no extrapolation needed by the combination scripts.
AEP_GRID = [0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01,
            0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001]

N_BOOTSTRAP = 5000
BOOT_SEED = 20260821
CONF_LOWER_PCT = 5
CONF_UPPER_PCT = 95

C_FIT = "#1a4f8a"
C_BAND = "#8fbcdb"
C_SYS = "#1a4f8a"
C_HIST = "#b7410e"

# ----------------------------------------------------------------------------


def read_systematic_usgs(path):
    """USGS peak-flow rdb -> {water_year: peak_cfs}, all rows in the file."""
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df.iloc[1:]                        # drop the "5s 15s 10d ..." spec row
    df = df[["peak_dt", "peak_va"]].copy()
    df["peak_dt"] = pd.to_datetime(df["peak_dt"], errors="coerce")
    df["peak_va"] = pd.to_numeric(df["peak_va"], errors="coerce")
    df = df.dropna(subset=["peak_dt", "peak_va"])
    # USGS water year: Oct-Dec belongs to the FOLLOWING calendar year's WY.
    df["water_year"] = df["peak_dt"].dt.year + (df["peak_dt"].dt.month >= 10).astype(int)
    return {int(r.water_year): float(r.peak_va) for r in df.itertuples()}


# A complete water year of 15-minute data is 365*96 = 35,040 records (more
# in a leap year). Below this, a file is a partial-year stub rather than a
# real water year -- caught once already: the 2006 file has 2,177 records
# (the Ecology gage started 08 Sep 2006, so "2006" is a ~3-week fragment,
# not a real annual peak) and the 2020 file cuts off after 4,848 records
# (~49 days). Both produced spuriously LOW "annual maxima" that blew up the
# skew estimate (-4.0 station skew on the first, unfiltered run) rather than
# erroring, which is exactly why this is a count check and not a try/except.
MIN_RECORDS_PER_WY = 30000       # ~83% coverage; tolerates real data gaps


def read_ecology_annual_max(path):
    """One Ecology 15-min FM file -> (water_year, peak_cfs), or None.

    Two header layouts are on disk (#Coweeman_Timing.py already documents
    this): a plain one starting straight at the column header, and a longer
    one with instructional text first. Both share the same column shape
    once the dashed divider line under the header is found: DATE, TIME,
    value, QUALITY, whitespace-separated.
    """
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("----"):
            start = i + 1
            break
    if start is None:
        return None

    best_val, best_date, n_records = None, None, 0
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        date_txt, val_txt = parts[0], parts[2]
        try:
            val = float(val_txt)
        except ValueError:
            continue
        n_records += 1
        if best_val is None or val > best_val:
            best_val, best_date = val, date_txt

    if best_val is None or n_records < MIN_RECORDS_PER_WY:
        return None

    month, year = int(best_date.split("/")[0]), int(best_date.split("/")[-1])
    water_year = year + 1 if month >= 10 else year
    return water_year, best_val


def read_systematic_ecology(pattern):
    out = {}
    for path in sorted(glob.glob(pattern)):
        result = read_ecology_annual_max(path)
        if result is None:
            continue
        wy, val = result
        # A file can straddle a water-year boundary if short; keep the
        # actual max found per water year rather than trusting the
        # filename's own year label.
        out[wy] = max(val, out.get(wy, 0.0))
    return out


def build_systematic_record():
    usgs = read_systematic_usgs(PEAKS_RDB)
    # 1996 is flagged peak_cd=7 in the rdb (a measured peak during an
    # otherwise-inactive stretch, not a normal water year of continuous
    # record) -- carried as a historic point instead, matching how CDID3
    # itself treats it.
    usgs = {wy: v for wy, v in usgs.items() if wy not in HISTORIC_PEAKS}
    ecology = read_systematic_ecology(ECOLOGY_GLOB)
    combined = dict(usgs)
    combined.update(ecology)     # the two periods do not overlap
    years = sorted(combined)
    values = np.array([combined[y] for y in years])
    return years, values


def fit_lp3(values):
    """Method-of-moments LP3 in log10 space: mean, std, station skew."""
    logs = np.log10(values)
    return logs.mean(), logs.std(ddof=1), stats.skew(logs, bias=False)


def station_skew_mse(skew, n):
    """Bulletin 17B (1982) MSE of a station skew estimate.

    MSE = 10 ** (A - B * log10(N/10))
      A = -0.33 + 0.08*|G|   if |G| <= 0.90   else  -0.52 + 0.30*|G|
      B =  0.94 - 0.26*|G|   if |G| <= 1.50   else   0.55
    """
    g = abs(skew)
    A = -0.33 + 0.08 * g if g <= 0.90 else -0.52 + 0.30 * g
    B = 0.94 - 0.26 * g if g <= 1.50 else 0.55
    return 10 ** (A - B * np.log10(n / 10.0))


def weighted_skew(station_skew, n):
    """Bulletin 17B weighted (regional-blended) skew, closed form."""
    mse_station = station_skew_mse(station_skew, n)
    return ((REGIONAL_SKEW_MSE * station_skew + mse_station * REGIONAL_SKEW)
             / (REGIONAL_SKEW_MSE + mse_station))


def lp3_quantiles(mean, std, skew, aep_grid):
    """AEP -> flow (cfs), Log-Pearson III via scipy's pearson3 frequency
    factor (falls back to the normal quantile at skew=0, which is what LP3
    reduces to anyway)."""
    aep = np.asarray(aep_grid, dtype=float)
    if abs(skew) < 1e-6:
        k = stats.norm.ppf(1 - aep)
    else:
        k = stats.pearson3.ppf(1 - aep, skew)
    return 10 ** (mean + k * std)


def bootstrap_band(values, aep_grid, n_boot=N_BOOTSTRAP, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    n = len(values)
    draws = np.empty((n_boot, len(aep_grid)))
    for i in range(n_boot):
        sample = rng.choice(values, size=n, replace=True)
        mean, std, skew = fit_lp3(sample)
        wskew = weighted_skew(skew, n)
        draws[i, :] = lp3_quantiles(mean, std, wskew, aep_grid)
    lower = np.percentile(draws, CONF_LOWER_PCT, axis=0)
    upper = np.percentile(draws, CONF_UPPER_PCT, axis=0)
    return lower, upper


def plotting_positions(values):
    """Weibull plotting position for the systematic record: rank/(n+1)."""
    n = len(values)
    order = np.argsort(values)[::-1]
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    return ranks / (n + 1.0)


def historic_plotting_positions(historic_dict, period_start, period_end):
    """Historic-adjusted plotting position for the two off-record points.

    Standard convention: rank the historic peaks among themselves and plot
    at rank/(H+1), where H is the length of the assumed historic period
    (from HISTORIC_PERIOD_START through the end of the systematic record),
    not the length of the systematic record itself. Plotting convention
    only -- does not feed the moment fit.
    """
    h = period_end - period_start + 1
    ordered = sorted(historic_dict.items(), key=lambda kv: -kv[1])
    return {yr: (rank + 1) / (h + 1.0) for rank, (yr, val) in enumerate(ordered)}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    years, values = build_systematic_record()
    n = len(values)
    print("Systematic Coweeman annual peaks: n=%d, %d-%d (gap 1985-2006 "
          "except the 1996 historic point)" % (n, min(years), max(years)))

    mean, std, station_skew = fit_lp3(values)
    wskew = weighted_skew(station_skew, n)
    print("Station skew=%.3f  regional-blended skew=%.3f" % (station_skew, wskew))

    best = lp3_quantiles(mean, std, wskew, AEP_GRID)
    lower, upper = bootstrap_band(values, AEP_GRID)

    out = pd.DataFrame({
        "AEP": AEP_GRID,
        "Value": best,
        "LowerConf": lower,
        "UpperConf": upper,
    })
    out["N_systematic"] = n
    out["StationSkew"] = station_skew
    out["WeightedSkew"] = wskew
    out["Method"] = "LP3_MoM_bootstrap_INTERIM_not_SSP"
    out.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV)
    print(out[["AEP", "Value", "LowerConf", "UpperConf"]].to_string(index=False))

    # -- plot --
    aep_pp = plotting_positions(values)
    hist_pp = historic_plotting_positions(HISTORIC_PEAKS, HISTORIC_PERIOD_START, max(years))

    fig, ax = plt.subplots(figsize=(9, 6.5))
    z_grid = stats.norm.ppf(1 - np.array(AEP_GRID))
    ax.plot(z_grid, best, color=C_FIT, lw=2, label="LP3 fit (method of moments, interim)")
    ax.fill_between(z_grid, lower, upper, color=C_BAND, alpha=0.35,
                     label="%d-%d%% bootstrap band" % (CONF_LOWER_PCT, CONF_UPPER_PCT))
    ax.scatter(stats.norm.ppf(1 - aep_pp), values, color=C_SYS, s=28, zorder=5,
               label="Systematic peaks (n=%d)" % n)
    for yr, p in hist_pp.items():
        ax.scatter(stats.norm.ppf(1 - p), HISTORIC_PEAKS[yr], color=C_HIST,
                   marker="^", s=70, zorder=6)
        ax.annotate(str(yr), (stats.norm.ppf(1 - p), HISTORIC_PEAKS[yr]),
                    textcoords="offset points", xytext=(6, 6), fontsize=9, color=C_HIST)
    ax.scatter([], [], color=C_HIST, marker="^", label="Historic peaks (shown, not fit)")

    ax.set_yscale("log")
    ax.set_xlabel("Standard normal variate  (z = Φ⁻¹(1 − AEP))")
    ax.set_ylabel("Coweeman River peak flow (cfs)")
    ax.set_title("Coweeman River near mouth — interim unregulated peak flow-frequency\n"
                  "(method-of-moments LP3 — NOT an HEC-SSP/EMA result, see script docstring)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    aep_ticks = [0.99, 0.5, 0.1, 0.02, 0.01, 0.002, 0.001, 0.0002, 0.0001]
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(stats.norm.ppf(1 - np.array(aep_ticks)))
    ax2.set_xticklabels(["%.2f%%" % (a * 100) for a in aep_ticks], rotation=45, fontsize=8)
    ax2.set_xlabel("AEP")

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("Wrote", PLOT_PNG)


if __name__ == "__main__":
    main()
