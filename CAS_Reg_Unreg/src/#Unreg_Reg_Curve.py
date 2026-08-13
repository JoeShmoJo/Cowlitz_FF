#Unreg_Reg_Curve.py
# -*- coding: utf-8 -*-
"""
Unregulated-regulated scatter and the inferred regulated frequency curve.

PART 1 -- SCATTER
    Unregulated peak against adjusted regulated peak at Castle Rock, drawn twice
    (arithmetic and log-log), both with a 1:1 line and the largest events called
    out by water year. The 1:1 line is the no-reservoir reference: points below
    it are the reduction the project achieved.

    The historic WCM_RC simulated pairs are drawn alongside, in their own
    colour. See "THE WCM_RC POINTS" below -- they are shown, not fitted.

PART 2 -- INFERRED REGULATED FREQUENCY CURVE
    The regulated curve is NOT fitted analytically. Regulated peaks do not
    follow an analytical distribution -- operating rules put hard breaks in
    them -- so the AEP is inherited from the unregulated side instead:

        unregulated AEP  ->  unregulated peak (Expected curve)
                         ->  unreg-reg relationship
                         ->  regulated peak at that AEP

    Peak-to-peak is used because peak and 1-day tied as the critical duration
    (log r-squared 0.850 vs 0.836, not distinguishable at n=44), and a
    peak-to-peak transform is the easier one to explain.

THE TRANSFORM IS A CENTRE-OF-MASS LINE, NOT A STRAIGHT LINE
------------------------------------------------------------
A single power law is a straight line in log-log space, and a straight line is
the wrong shape here. Regulated peaks are shaped by operating rules: the
reservoir holds the small events back hard, loses ground through the middle,
and at the top the relationship flattens towards the unregulated flow as the
project runs out of storage. Forcing one slope through all of that puts the
line above the data in one range and below it in another, and the error lands
where it matters most -- the upper end.

TRANSFORM_METHOD = "loess" therefore draws a locally weighted centre-of-mass
line through the scatter instead: at every point on the curve the slope comes
from the nearby data only (tricube weights, LOESS_SPAN of the sample), so the
line is free to bend where the data bends. It is forced to increase
monotonically, and optionally clipped at the 1:1 line, since a regulated peak
above its own unregulated peak is not physical.

TRANSFORM_METHOD = "power" restores the old single power law. It is still drawn
as a thin reference line either way, so the difference is visible on the plot.

Both are rough lines through a scatter, meant to be replaced by a hand-drawn
curve once the synthetics populate the upper end.

THE WCM_RC POINTS
-----------------
The WCM_RC ResSim run simulates BOTH the regulated and the unregulated Castle
Rock flow over the historic period, so it supplies a simulated pair for every
water year -- including the years the adjusted record screens out, and the
pre-regulation years the adjusted record cannot reach at all.

They are drawn in their own colour and, by default, are NOT fitted, because
they cannot be corrected the way the adjusted record can: the adjustment is
built from the difference between the WCM_RC and Obs_RC runs, and where there
is no Obs_RC data there is no correction to apply. They are simulated flows
carrying whatever bias the model carries, with nothing observed to anchor them.
Treat them as context for the shape of the relationship, not as evidence for
its position. Set INCLUDE_WCM_IN_FIT = True to fit them anyway -- and say so in
the memo if you do.

The same reg-above-unreg screen used in #Adjusted_Peak_Record.py is applied to
them, at the same threshold, so a non-physical simulated pair is not shown as
though it were usable.

PART 3 -- COMPARISON WITH THE 2009 STUDY
    The adopted regulated frequency curve from the 2009 study is hard coded
    below and drawn on the same axes, with the differences tabulated in the
    output CSV. It is the only external check available on the upper end.

CAUTION
    The transform is supported over the observed unregulated range only. Beyond
    that it is extrapolation whichever method is used, and the plot marks where
    the support stops.
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
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
DATASET_CSV = r"../output/diagnostics/critical_duration_adjusted_dataset.csv"
FITS_CSV = r"../output/critical_duration_adjusted_fits.csv"
UNREG_FREQ_CSV = r"../../CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv"
# Screening authority. Any water year not marked eligible here is dropped, even
# if it is present in DATASET_CSV -- that file can be stale.
ADJUSTED_PEAKS_CSV = r"../output/adjusted_peaks.csv"
ENFORCE_SCREENING = True

# Historic simulated pairs from the WCM_RC run, written by
# #Extract_Ensemble_To_Timeseries.py with SET_NAME = "ResSim_WCM_RC".
WCM_WY_CSV = r"../output/diagnostics/ResSim_WCM_RC_reg_vs_unreg_wy.csv"
WCM_PART_B = "CastleRock_NWS"
SHOW_WCM_POINTS = True
# Fit them as well? Default False -- there is no Obs_RC data behind them, so
# there is no adjustment available to correct them.
INCLUDE_WCM_IN_FIT = False
# Restrict them to the regulated era (e.g. 1974) or None for the whole record.
WCM_FIRST_WY = None
# Apply the same reg-above-unreg screen used on the adjusted record.
WCM_SCREEN_REG_OVER_UNREG = True
WCM_REG_OVER_UNREG_THRESHOLD_CFS = 60000.0

OUT_DIR = r"../output"
PLOT_STEM = r"../output/diagnostics/unreg_reg"

UNREG_COL = "unreg_Peak_1-hr"      # peak-to-peak transform
FREQ_DURATION = "Peak"
FREQ_VALUE_COL = "Expected"        # the expected-probability curve

CALLOUT_TOP_N = 8                  # label this many largest events
CALLOUT_MIN_CFS = 90000.0          # ...and anything above this

# --- how the unreg -> reg transform is drawn --------------------------------
#   "loess" : locally weighted centre-of-mass line through the scatter, free to
#             bend. The default, and the reason is in the docstring.
#   "power" : single power law reg = a * unreg^b, a straight line in log-log.
TRANSFORM_METHOD = "loess"
# Fraction of the sample inside each local window. Larger = smoother and
# straighter; smaller = follows the data more closely and gets noisier. 0.5-0.8
# is a sensible range at n = 40-ish.
LOESS_SPAN = 0.65
# Never let the drawn curve decrease as the unregulated flow increases.
ENFORCE_MONOTONIC = True
# Never let the drawn regulated flow exceed the unregulated flow it came from.
CLIP_TO_UNREG = True
# Draw the single power law as a thin reference line for comparison.
SHOW_POWER_LAW_REFERENCE = True

# --- 2009 study adopted regulated frequency curve ----------------------------
# AEP in percent, discharge in cfs. The external check on the upper end.
CURVE_2009_LABEL = "2009 study, regulated"
CURVE_2009 = [
    (99.0, 21700.0), (95.0, 28100.0), (90.0, 32400.0), (80.0, 38800.0),
    (70.0, 44400.0), (60.0, 49800.0), (50.0, 55500.0), (40.0, 59500.0),
    (30.0, 64100.0), (20.0, 70000.0), (10.0, 74000.0), (5.0, 79000.0),
    (4.0, 81200.0), (2.0, 88000.0), (1.0, 97000.0), (0.7, 104000.0),
    (0.5, 110000.0), (0.2, 156000.0), (0.1, 240000.0), (0.08, 270000.0),
    (0.05, 300000.0), (0.01, 390000.0),
]

# How the observed points are placed on the frequency plot.
#   "from_curve" : each year's AEP is read off the UNREGULATED frequency curve
#                  at its own unregulated peak, and the regulated peak is drawn
#                  at that same AEP. This is the method itself made visible --
#                  the regulated point inherits the unregulated AEP -- and it is
#                  the only consistent choice here, because these years are the
#                  regulated-era subset of a 98-year record. Weibull positions
#                  computed on the subset would place them at far higher AEPs
#                  than the curve, which is fitted to all 98.
#   "weibull"    : plotting positions from this subset alone. Kept for
#                  comparison; expect the points to sit left of the curve.
PLOTTING_BASIS = "from_curve"

# SSP-style frequency axis
AEP_TICKS = [0.999, 0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02,
             0.01, 0.005, 0.002, 0.001]
AEP_LIMITS = (0.999, 0.0005)
FLOW_LIMITS = (10000.0, 400000.0)

# Colours kept in one place so the scatter and the frequency plot agree.
C_UNREG = "#2c7fb8"
C_REG = "#c0392b"
C_WCM = "#7d3c98"          # WCM_RC simulated pairs -- deliberately distinct
C_2009 = "#117a65"

# ----------------------------------------------------------------------------


def weibull_plotting_positions(values):
    """Weibull plotting positions (i/(n+1)) as AEP, largest first."""
    clean = np.sort(np.asarray(values, dtype=float))[::-1]
    n = len(clean)
    return clean, (np.arange(1, n + 1)) / float(n + 1)


def fit_power_law(x, y):
    """log10(y) on log10(x). Returns a, b, r2 and the log-space standard error."""
    good = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    lx, ly = np.log10(x[good]), np.log10(y[good])
    fit = stats.linregress(lx, ly)
    resid = ly - (fit.slope * lx + fit.intercept)
    return {"a": 10 ** fit.intercept, "b": fit.slope, "r2": fit.rvalue ** 2,
            "se_dex": float(np.std(resid, ddof=2)), "n": int(good.sum()),
            "x_min": float(x[good].min()), "x_max": float(x[good].max())}


def apply_power_law(fit, x):
    return fit["a"] * np.asarray(x, dtype=float) ** fit["b"]


def loess_at(lx, ly, x0, span):
    """One locally weighted linear estimate at x0, all in log10 space.

    Tricube weights over the nearest ceil(span*n) points. Local LINEAR rather
    than local mean, so the line keeps a sensible slope at the ends of the data
    and where it is extrapolating.
    """
    n = len(lx)
    k = int(np.ceil(span * n))
    k = max(4, min(k, n))
    dist = np.abs(lx - x0)
    near = np.argsort(dist)[:k]
    d = dist[near]
    d_max = d.max()
    weights = np.ones(k) if d_max <= 0 else (1.0 - (d / d_max) ** 3) ** 3
    weights = np.clip(weights, 1e-8, None)
    xs, ys = lx[near], ly[near]
    w_sum = weights.sum()
    x_bar = (weights * xs).sum() / w_sum
    y_bar = (weights * ys).sum() / w_sum
    sxx = (weights * (xs - x_bar) ** 2).sum()
    if sxx <= 0:
        return y_bar
    slope = (weights * (xs - x_bar) * (ys - y_bar)).sum() / sxx
    return y_bar + slope * (x0 - x_bar)


def build_transform(x, y, method):
    """The unreg -> reg relationship. Returns a dict usable by apply_transform.

    Always carries the single power law too, so it can be drawn as a reference
    and reported in the log even when the LOESS line is the one adopted.
    """
    good = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[good], np.asarray(y)[good]
    power = fit_power_law(x, y)
    fit = {"method": method, "power": power, "n": int(good.sum()),
           "x_min": float(x.min()), "x_max": float(x.max()),
           "x_obs": x, "y_obs": y}
    if method == "power":
        fit["se_dex"] = power["se_dex"]
        fit["r2"] = power["r2"]
        return fit

    lx, ly = np.log10(x), np.log10(y)
    fit["lx"], fit["ly"] = lx, ly
    fit["span"] = LOESS_SPAN
    predicted = np.array([loess_at(lx, ly, v, LOESS_SPAN) for v in lx])
    resid = ly - predicted
    # Effective parameters of a LOESS fit are not 2; ddof=2 is a rough and
    # slightly optimistic stand-in, matching what the power law reports so the
    # two bands are comparable.
    fit["se_dex"] = float(np.std(resid, ddof=2))
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    fit["r2"] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return fit


def apply_transform(fit, x):
    """Regulated flow implied by an unregulated flow."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if fit["method"] == "power":
        out = apply_power_law(fit["power"], x)
    else:
        lx = np.log10(np.clip(x, 1e-6, None))
        out = np.array([10.0 ** loess_at(fit["lx"], fit["ly"], v, fit["span"])
                        for v in lx])
    if ENFORCE_MONOTONIC and len(out) > 1:
        order = np.argsort(x)
        ordered = np.maximum.accumulate(out[order])
        result = np.empty_like(out)
        result[order] = ordered
        out = result
    if CLIP_TO_UNREG:
        out = np.minimum(out, x)
    return out


def transform_label(fit):
    if fit["method"] == "power":
        p = fit["power"]
        return ("Power law: reg = %.4g x unreg$^{%.3f}$  (r$^2$=%.3f)"
                % (p["a"], p["b"], p["r2"]))
    return ("LOESS centre of mass (span %.2f, r$^2$=%.3f)"
            % (fit["span"], fit["r2"]))


def load_eligible_wys(csv_path):
    """Water years that passed screening in #Adjusted_Peak_Record.py."""
    if not (ENFORCE_SCREENING and os.path.isfile(csv_path)):
        return None
    table = pd.read_csv(csv_path)
    if "screen_passed" not in table.columns:
        print("WARNING %s has no screen_passed column -- screening NOT enforced."
              % csv_path)
        print("        Re-run #Adjusted_Peak_Record.py; screened-out years may")
        print("        be sitting in the scatter and the frequency plot.")
        return None
    return table


def load_wcm_points(csv_path):
    """Historic simulated (unreg, reg) pairs from the WCM_RC run.

    Returns (kept, dropped). Nothing is fitted from these by default; see the
    docstring.
    """
    if not (SHOW_WCM_POINTS and os.path.isfile(csv_path)):
        if SHOW_WCM_POINTS:
            print("WCM_RC  not found: %s -- simulated points not drawn" % csv_path)
        return None, None
    table = pd.read_csv(csv_path)
    if "part_b" in table.columns:
        table = table[table["part_b"] == WCM_PART_B]
    table = table.dropna(subset=["reg_peak", "unreg_peak"]).copy()
    if WCM_FIRST_WY is not None:
        table = table[table["WY"] >= WCM_FIRST_WY]
    dropped = table.iloc[0:0]
    if WCM_SCREEN_REG_OVER_UNREG:
        over = ((table["reg_peak"] - table["unreg_peak"]) > 1.0) & \
               (table["reg_peak"] >= WCM_REG_OVER_UNREG_THRESHOLD_CFS)
        dropped = table[over]
        table = table[~over]
    return table.reset_index(drop=True), dropped.reset_index(drop=True)


def curve_2009_frame():
    """The 2009 adopted regulated curve as AEP (fraction) and cfs."""
    table = pd.DataFrame(CURVE_2009, columns=["AEP_pct", "cfs"])
    table["AEP"] = table["AEP_pct"] / 100.0
    return table.sort_values("AEP", ascending=False).reset_index(drop=True)


def interp_2009(aep, table_2009):
    """2009 discharge at any AEP: log10(Q) interpolated against the z-variate.

    Interpolating in z and log-flow is the same space the curve is drawn in, so
    the interpolated values sit on the plotted line rather than cutting corners
    off it. Outside the tabulated range the result is NaN, not an extrapolation.
    """
    z_known = stats.norm.ppf(1.0 - table_2009["AEP"].values)
    q_known = np.log10(table_2009["cfs"].values)
    order = np.argsort(z_known)
    z = stats.norm.ppf(1.0 - np.asarray(aep, dtype=float))
    out = np.interp(z, z_known[order], q_known[order], left=np.nan, right=np.nan)
    return 10.0 ** out


def annotate_points(ax, x, y, labels):
    """Callouts placed alternately above and below to reduce collisions."""
    for i, (xi, yi, text) in enumerate(zip(x, y, labels)):
        dy = 14 if i % 2 == 0 else -18
        ax.annotate(text, (xi, yi), xytext=(10, dy), textcoords="offset points",
                    fontsize=8, color=C_REG,
                    arrowprops=dict(arrowstyle="-", color=C_REG, lw=0.7))


def plot_scatter(data, fit, wcm, stem):
    """Arithmetic and log-log scatter, each with a 1:1 line and callouts."""
    x = data[UNREG_COL].values
    y = data["reg_peak"].values
    big = data.nlargest(CALLOUT_TOP_N, UNREG_COL)
    big = pd.concat([big, data[data[UNREG_COL] >= CALLOUT_MIN_CFS]]).drop_duplicates()

    for log_axes in (False, True):
        fig, ax = plt.subplots(figsize=(9.5, 8.5))

        lo = [x.min(), y.min()]
        hi = [x.max(), y.max()]
        if wcm is not None and len(wcm):
            ax.scatter(wcm["unreg_peak"], wcm["reg_peak"], s=30, marker="D",
                       facecolor="none", edgecolor=C_WCM, lw=0.9, zorder=2,
                       label="WCM_RC simulated pair (no Obs_RC correction, "
                             "not fitted)" if not INCLUDE_WCM_IN_FIT
                             else "WCM_RC simulated pair (fitted)")
            lo += [wcm["unreg_peak"].min(), wcm["reg_peak"].min()]
            hi += [wcm["unreg_peak"].max(), wcm["reg_peak"].max()]

        ax.scatter(x, y, s=46, facecolor=C_UNREG, edgecolor="0.25", lw=0.6,
                   zorder=3, label="Water year (adjusted regulated peak)")
        lim = [min(lo) * 0.85, max(hi) * 1.15]
        ax.plot(lim, lim, color="k", lw=1.2, ls="--", zorder=2,
                label="1:1 (no reservoir effect)")

        xs = np.geomspace(max(lim[0], 1.0), lim[1], 300)
        centre = apply_transform(fit, xs)
        ax.plot(xs, centre, color=C_REG, lw=1.8, zorder=4,
                label=transform_label(fit))
        band = 10 ** fit["se_dex"]
        ax.fill_between(xs, centre / band, centre * band, color=C_REG,
                        alpha=0.12, zorder=1, label="+/- 1 std error")
        if SHOW_POWER_LAW_REFERENCE and fit["method"] != "power":
            ax.plot(xs, apply_power_law(fit["power"], xs), color=C_REG, lw=1.0,
                    ls=":", zorder=3,
                    label="Single power law, for reference (r$^2$=%.3f)"
                          % fit["power"]["r2"])

        annotate_points(ax, big[UNREG_COL].values, big["reg_peak"].values,
                        ["WY%d" % w for w in big["WY"]])
        if log_axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.grid(which="both", alpha=0.25)
        else:
            ax.set_xlim(0, lim[1])
            ax.set_ylim(0, lim[1])
            ax.grid(alpha=0.3)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
        ax.set_xlabel("Unregulated peak at Castle Rock (cfs)")
        ax.set_ylabel("Regulated peak at Castle Rock (cfs)")
        ax.set_title("Castle Rock unregulated vs regulated peak%s\n"
                     "points below the 1:1 line are the reduction the project achieved"
                     % ("  (log-log)" if log_axes else ""), fontsize=11)
        ax.legend(loc="upper left", fontsize=8.5)
        fig.tight_layout()
        fig.savefig("%s_scatter_%s.png" % (stem, "loglog" if log_axes else "linear"),
                    dpi=150)
        plt.close(fig)


def probability_axis(ax, ticks, limits):
    """SSP-style normal-probability axis, AEP decreasing to the right."""
    z_ticks = stats.norm.ppf(1.0 - np.array(ticks))
    ax.set_xlim(stats.norm.ppf(1.0 - limits[0]), stats.norm.ppf(1.0 - limits[1]))
    ax.xaxis.set_major_locator(FixedLocator(z_ticks))
    ax.set_xticklabels(["%g" % (t * 100) for t in ticks])
    ax.set_xlabel("Annual exceedance probability (%)")


def aep_from_unreg_curve(values, unreg_curve, aep_values):
    """AEP of a flow, read off the unregulated frequency curve."""
    order = np.argsort(unreg_curve)
    return np.interp(np.asarray(values, dtype=float), unreg_curve[order],
                     aep_values[order])


def plot_frequency(freq, data, fit, wcm, table_2009, stem):
    """Unregulated and inferred regulated frequency curves, SSP idiom."""
    fig, ax = plt.subplots(figsize=(11.5, 8.8))

    aep_values = freq["AEP"].values
    z = stats.norm.ppf(1.0 - aep_values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    reg_curve = apply_transform(fit, unreg_curve)
    band = 10 ** fit["se_dex"]

    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.2, zorder=4,
            label="Unregulated peak (%s curve)" % FREQ_VALUE_COL.lower())
    ax.plot(z, reg_curve, color=C_REG, lw=2.2, zorder=4,
            label="Regulated peak, inferred (%s)"
                  % ("LOESS centre of mass" if fit["method"] == "loess"
                     else "power law"))
    ax.fill_between(z, reg_curve / band, reg_curve * band, color=C_REG,
                    alpha=0.14, zorder=1, label="Regulated, +/- 1 std error")
    if SHOW_POWER_LAW_REFERENCE and fit["method"] != "power":
        ax.plot(z, apply_power_law(fit["power"], unreg_curve), color=C_REG,
                lw=1.0, ls=":", zorder=4,
                label="Regulated via the single power law (reference)")

    # --- the 2009 adopted regulated curve -----------------------------------
    if table_2009 is not None and len(table_2009):
        z_2009 = stats.norm.ppf(1.0 - table_2009["AEP"].values)
        ax.plot(z_2009, table_2009["cfs"].values, color=C_2009, lw=1.9,
                ls="--", marker="o", ms=4, mfc="white", zorder=5,
                label=CURVE_2009_LABEL)

    # --- observed / adjusted points -----------------------------------------
    if PLOTTING_BASIS == "from_curve":
        # AEP of each year's UNREGULATED peak, read off the unregulated curve;
        # the regulated peak of that same year is drawn at the same AEP.
        aep_of = aep_from_unreg_curve(data[UNREG_COL].values, unreg_curve,
                                      aep_values)
        zz = stats.norm.ppf(1.0 - np.clip(aep_of, 1e-6, 1 - 1e-6))
        ax.plot(zz, data[UNREG_COL].values, ls="none", marker="o", ms=5,
                mfc=C_UNREG, mec="0.2", mew=0.5, zorder=6,
                label="Unregulated peaks (AEP from the curve)")
        ax.plot(zz, data["reg_peak"].values, ls="none", marker="s", ms=5,
                mfc=C_REG, mec="0.2", mew=0.5, zorder=6,
                label="Adjusted regulated peaks (same AEP as their unreg peak)")
        for i in range(len(zz)):
            ax.plot([zz[i], zz[i]],
                    [data[UNREG_COL].values[i], data["reg_peak"].values[i]],
                    color="0.6", lw=0.5, zorder=3)
        # --- WCM_RC simulated pairs, same treatment, own colour -------------
        if wcm is not None and len(wcm):
            aep_w = aep_from_unreg_curve(wcm["unreg_peak"].values, unreg_curve,
                                         aep_values)
            zw = stats.norm.ppf(1.0 - np.clip(aep_w, 1e-6, 1 - 1e-6))
            ax.plot(zw, wcm["reg_peak"].values, ls="none", marker="D", ms=4.5,
                    mfc="none", mec=C_WCM, mew=1.0, zorder=5,
                    label="WCM_RC simulated regulated peaks "
                          "(no Obs_RC correction)")
    else:
        uv, ua = weibull_plotting_positions(data[UNREG_COL].values)
        rv, ra = weibull_plotting_positions(data["reg_peak"].values)
        ax.plot(stats.norm.ppf(1.0 - ua), uv, ls="none", marker="o", ms=5,
                mfc=C_UNREG, mec="0.2", mew=0.5, zorder=6,
                label="Unregulated peaks (Weibull, n=%d subset)" % len(uv))
        ax.plot(stats.norm.ppf(1.0 - ra), rv, ls="none", marker="s", ms=5,
                mfc=C_REG, mec="0.2", mew=0.5, zorder=6,
                label="Adjusted regulated peaks (Weibull, n=%d subset)" % len(rv))
        if wcm is not None and len(wcm):
            wv, wa = weibull_plotting_positions(wcm["reg_peak"].values)
            ax.plot(stats.norm.ppf(1.0 - wa), wv, ls="none", marker="D", ms=4.5,
                    mfc="none", mec=C_WCM, mew=1.0, zorder=5,
                    label="WCM_RC simulated regulated peaks (Weibull, n=%d)"
                          % len(wv))

    # where the transform stops being supported by data
    supported = unreg_curve <= fit["x_max"]
    if supported.any() and not supported.all():
        z_edge = z[supported].max()
        ax.axvline(z_edge, color="0.35", lw=1.0, ls=":", zorder=2)
        ax.text(z_edge, FLOW_LIMITS[1] * 0.92,
                "  transform supported to here\n  (unreg %s cfs)"
                % format(int(fit["x_max"]), ","),
                fontsize=8, color="0.3", va="top")

    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    probability_axis(ax, AEP_TICKS, AEP_LIMITS)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("Castle Rock peak flow frequency\n"
                 "regulated AEP inherited from the unregulated curve through a "
                 "%s relationship"
                 % ("centre-of-mass (LOESS)" if fit["method"] == "loess"
                    else "power-law"), fontsize=12)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92)
    if wcm is not None and len(wcm):
        ax.text(0.995, 0.015,
                "WCM_RC points are simulated on both axes and carry no Obs_RC\n"
                "correction -- shown for shape, not used to place the curve.",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color=C_WCM,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=C_WCM,
                          alpha=0.85, lw=0.8))
    fig.tight_layout()
    fig.savefig("%s_frequency.png" % stem, dpi=150)
    plt.close(fig)
    return reg_curve


def plot_2009_comparison(out, stem):
    """Ratio of the inferred regulated curve to the 2009 adopted curve."""
    sub = out.dropna(subset=["reg_2009_cfs"])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(11, 5.2))
    z = stats.norm.ppf(1.0 - sub["AEP"].values)
    ax.plot(z, 100.0 * (sub["reg_inferred_cfs"] / sub["reg_2009_cfs"] - 1.0),
            color=C_REG, lw=2.0, marker="o", ms=3.5)
    ax.axhline(0.0, color=C_2009, lw=1.4, ls="--")
    ax.text(z.min(), 1.0, " 2009 adopted curve", color=C_2009, fontsize=9,
            va="bottom")
    probability_axis(ax, AEP_TICKS, AEP_LIMITS)
    ax.set_ylabel("2026 inferred vs 2009 (%)")
    ax.set_title("Inferred regulated curve against the 2009 study\n"
                 "positive = the 2026 estimate is higher", fontsize=11)
    ax.grid(alpha=0.35)
    fig.tight_layout()
    fig.savefig("%s_2009_comparison.png" % stem, dpi=150)
    plt.close(fig)


def main():
    for path in (OUT_DIR, os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    data = pd.read_csv(DATASET_CSV).dropna(subset=[UNREG_COL, "reg_peak"])

    # --- screening is enforced here as well -----------------------------------
    # DATASET_CSV is written by #Critical_Duration_Adjusted.py, which already
    # filters on screen_passed -- but only as of the run that produced it. A
    # stale dataset is the way a screened year gets back into these plots, so
    # the authority is re-read and applied rather than trusted.
    screening = load_eligible_wys(ADJUSTED_PEAKS_CSV)
    if screening is not None:
        eligible = set(screening.loc[screening["screen_passed"].astype(bool), "WY"])
        removed = sorted(set(data["WY"]) - eligible)
        if removed:
            codes = screening.set_index("WY")
            print("Screening : dropped %d water year(s) from %s -- %s"
                  % (len(removed), os.path.basename(DATASET_CSV),
                     ", ".join("WY%d (%s)"
                               % (w, codes.loc[w, "screen_code"]
                                  if "screen_code" in codes.columns else "screened")
                               for w in removed)))
            print("            re-run #Critical_Duration_Adjusted.py to refresh it")
        data = data[data["WY"].isin(eligible)].reset_index(drop=True)

    wcm, wcm_dropped = load_wcm_points(WCM_WY_CSV)
    if wcm is not None:
        print("WCM_RC    : %d simulated pairs (WY%d..%d)%s"
              % (len(wcm), int(wcm["WY"].min()), int(wcm["WY"].max()),
                 ", %d screened out (reg above unreg at >= %s cfs)"
                 % (len(wcm_dropped), format(int(WCM_REG_OVER_UNREG_THRESHOLD_CFS), ","))
                 if wcm_dropped is not None and len(wcm_dropped) else ""))
        if wcm_dropped is not None and len(wcm_dropped):
            print("            screened: %s"
                  % ", ".join("WY%d" % w for w in wcm_dropped["WY"]))

    fit_x = data[UNREG_COL].values
    fit_y = data["reg_peak"].values
    if INCLUDE_WCM_IN_FIT and wcm is not None and len(wcm):
        fit_x = np.concatenate([fit_x, wcm["unreg_peak"].values])
        fit_y = np.concatenate([fit_y, wcm["reg_peak"].values])
    fit = build_transform(fit_x, fit_y, TRANSFORM_METHOD)

    freq = pd.read_csv(UNREG_FREQ_CSV)
    freq = freq[freq["Duration"] == FREQ_DURATION].sort_values("AEP", ascending=False)
    freq = freq.dropna(subset=["AEP", FREQ_VALUE_COL])
    table_2009 = curve_2009_frame()

    print("=" * 78)
    print("Transform : %s   n = %d%s"
          % (TRANSFORM_METHOD, fit["n"],
             "  (WCM_RC pairs INCLUDED in the fit)" if INCLUDE_WCM_IN_FIT else ""))
    if TRANSFORM_METHOD == "loess":
        print("            LOESS span %.2f, pseudo r2 = %.4f" % (LOESS_SPAN, fit["r2"]))
        print("            power law for reference: reg = %.4g x unreg^%.4f, r2 = %.4f"
              % (fit["power"]["a"], fit["power"]["b"], fit["power"]["r2"]))
    else:
        print("            reg = %.4g x unreg^%.4f   log r2 = %.4f"
              % (fit["power"]["a"], fit["power"]["b"], fit["power"]["r2"]))
    print("Supported over unregulated %s to %s cfs"
          % (format(int(fit["x_min"]), ","), format(int(fit["x_max"]), ",")))
    print("Scatter about the line: x/ %.3f (1 sigma)" % 10 ** fit["se_dex"])
    print("Monotonic enforced: %s   clipped at 1:1: %s"
          % (ENFORCE_MONOTONIC, CLIP_TO_UNREG))
    print("=" * 78)

    plot_scatter(data, fit, wcm, PLOT_STEM)
    reg_curve = plot_frequency(freq, data, fit, wcm, table_2009, PLOT_STEM)

    out = freq[["AEP", "Value", FREQ_VALUE_COL]].copy()
    out = out.rename(columns={"Value": "unreg_computed_cfs",
                              FREQ_VALUE_COL: "unreg_expected_cfs"})
    band = 10 ** fit["se_dex"]
    out["reg_inferred_cfs"] = reg_curve
    out["reg_lower_1se_cfs"] = reg_curve / band
    out["reg_upper_1se_cfs"] = reg_curve * band
    out["reg_powerlaw_cfs"] = apply_power_law(fit["power"],
                                              out["unreg_expected_cfs"].values)
    out["reduction_pct"] = 100.0 * (1.0 - reg_curve / out["unreg_expected_cfs"])
    out["extrapolated"] = out["unreg_expected_cfs"] > fit["x_max"]
    out["reg_2009_cfs"] = interp_2009(out["AEP"].values, table_2009)
    out["reg_minus_2009_cfs"] = out["reg_inferred_cfs"] - out["reg_2009_cfs"]
    out["reg_vs_2009_pct"] = 100.0 * (out["reg_inferred_cfs"] / out["reg_2009_cfs"]
                                      - 1.0)
    out.to_csv(os.path.join(OUT_DIR, "regulated_frequency_inferred.csv"),
               index=False, float_format="%.1f")

    plot_2009_comparison(out, PLOT_STEM)

    show = out[out["AEP"].isin([0.5, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002])]
    print("Observed points placed by: %s" % PLOTTING_BASIS)
    print("\nINFERRED REGULATED CURVE")
    print(show[["AEP", "unreg_expected_cfs", "reg_inferred_cfs",
                "reg_lower_1se_cfs", "reg_upper_1se_cfs", "reduction_pct",
                "extrapolated"]]
          .round({"AEP": 4, "unreg_expected_cfs": 0, "reg_inferred_cfs": 0,
                  "reg_lower_1se_cfs": 0, "reg_upper_1se_cfs": 0,
                  "reduction_pct": 1}).to_string(index=False))

    # --- comparison with the 2009 study --------------------------------------
    print("\nAGAINST THE 2009 STUDY (regulated)")
    at_2009 = out.dropna(subset=["reg_2009_cfs"])
    if at_2009.empty:
        print("   no overlap between the AEPs on the 2026 curve and the 2009 table")
    else:
        wanted = [0.99, 0.5, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002, 0.001]
        rows = at_2009[at_2009["AEP"].isin(wanted)]
        if rows.empty:
            rows = at_2009
        print(rows[["AEP", "reg_inferred_cfs", "reg_2009_cfs",
                    "reg_minus_2009_cfs", "reg_vs_2009_pct"]]
              .round({"AEP": 4, "reg_inferred_cfs": 0, "reg_2009_cfs": 0,
                      "reg_minus_2009_cfs": 0, "reg_vs_2009_pct": 1})
              .to_string(index=False))
        print("   median difference across the overlap: %+.1f%%   range %+.1f%% to %+.1f%%"
              % (at_2009["reg_vs_2009_pct"].median(),
                 at_2009["reg_vs_2009_pct"].min(),
                 at_2009["reg_vs_2009_pct"].max()))
        print("   The 2009 curve turns sharply upward above the 0.2% AEP "
              "(156,000 -> 390,000 cfs).")
        print("   Any comparison out there is between two extrapolations, not "
              "between two records.")

    n_ex = int(out["extrapolated"].sum())
    if n_ex:
        print("\n   %d of %d AEPs are beyond the supported range -- these depend on"
              % (n_ex, len(out)))
        print("   the relationship holding past any event on record. The synthetics")
        print("   exist to replace that extrapolation with simulated points.")
    print("-" * 78)
    print("Plots : %s_scatter_linear.png, %s_scatter_loglog.png,"
          % (PLOT_STEM, PLOT_STEM))
    print("        %s_frequency.png, %s_2009_comparison.png" % (PLOT_STEM, PLOT_STEM))
    print("Table : %s/regulated_frequency_inferred.csv" % OUT_DIR)


main()
