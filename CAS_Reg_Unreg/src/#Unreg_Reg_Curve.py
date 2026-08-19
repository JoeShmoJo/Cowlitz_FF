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
its position. FIT_SOURCES["wcm_rc"] = True fits them anyway -- and say so in the
memo if you do.

The same reg-above-unreg screen used in #Adjusted_Peak_Record.py is applied to
them, at the same threshold, so a non-physical simulated pair is not shown as
though it were usable.

PART 3 -- COMPARISON WITH THE 2009 STUDY
    The adopted regulated frequency curve from the 2009 study is hard coded
    below and drawn on the same axes, with the differences tabulated in the
    output CSV. It is the only external check available on the upper end.

    The 2009 curve is a DIAGNOSTIC comparison, not part of the result. The
    adopted report figure (FINAL_SHOW_2009 = False) leaves it off.

PART 4 -- UNCERTAINTY ON THE REGULATED CURVE
    Follows EM 1110-2-1619 (29 Sep 2025) Sec 4-4.b(3) eq 4-6, applied to flow
    instead of stage: independent uncertainty sources combine by adding their
    VARIANCES, i.e. S_total = sqrt(S_1^2 + ... + S_n^2). Two sources here:

      1. FLOW FREQUENCY. The unregulated quantile at an AEP is itself an
         estimate from a 98-year record. HEC-SSP reports this as confidence
         limits about the computed curve.
      2. THE UNREG-REG TRANSFORM. At one unregulated magnitude the regulated
         peak still depends on the shape of the flood, which is what the
         scatter of the (unreg, reg) pairs about the transform line measures.

    Each is reduced to its OWN standard deviation first (frequency_sigma_dex,
    transform_sigma_dex_split), in whatever units are natural for that source
    -- log10(cfs) for both, since that is the space HEC-SSP's own limits and
    the LOESS residuals are already in. The two sigmas are combined by eq 4-6,
    and only THEN is a single z (from UNCERTAINTY_CONF_LEVEL) applied to the
    combined sigma -- not one z per source. Combining "distance to a chosen
    percentile" directly, without ever naming sigma, is the same arithmetic
    (z factors out of a root-sum-of-squares), but only if every source's
    distance was measured at the SAME percentile; that is why this script does
    not do it, since the frequency term is natively reported at SSP's own 90%
    while the transform term is estimated from the data at whatever precision
    the local scatter supports.

    THE FREQUENCY TERM HAS TO BE IN REGULATED CFS, NOT UNREGULATED CFS, before
    it can be combined with the transform term -- summing an unregulated-flow
    delta with a regulated-flow delta is a dimensional error. This is done by
    literally pushing the unregulated flow AT ITS OWN UPPER/LOWER BOUND through
    the SAME best-fit transform used for the central curve, and taking the gap
    from the regulated best estimate:

        freq_term_hi = Transform(Unreg_upper) - Transform(Unreg_best)
        freq_term_lo = Transform(Unreg_best) - Transform(Unreg_lower)

    where Unreg_upper/lower are the unregulated flow at the chosen z, from its
    own asymmetric sigma (frequency_sigma_dex). This evaluates the curve's
    actual bend rather than a linear (constant-slope) approximation to it --
    simpler to explain, and more correct when the band is wide, since the
    transform bends noticeably between the best estimate and its confidence
    limit at the high end (see transform_log_slope, kept only for reporting
    context: the local slope runs about 0.57 through the middle, where the
    reservoir absorbs most of an increase, and about 1.5 at the top, where the
    transform bends toward pass-through -- which is why the frequency term
    below grows disproportionately large there).

    THE TRANSFORM TERM is independent of the one above: held at the
    unregulated BEST ESTIMATE (not its bound), how far do actual regulated
    points scatter above/below the fitted line there, expressed in regulated
    cfs the same way -- reg_best * (10**(z*sigma) - 1).

    COMBINED ADDITIVELY, IN CFS, PER SIDE -- this is eq 4-6 with S already
    converted to cfs at the chosen z, and it is the literal form the guidance
    was given back to this project in:

        Upper = RegulatedBest + sqrt(freq_term_hi^2 + transform_term_hi^2)
        Lower = RegulatedBest - sqrt(freq_term_lo^2 + transform_term_lo^2)

    computed independently for each side, which is how the asymmetry survives
    without ever needing to be a special case.

    THE ASYMMETRY, ON BOTH SOURCES. The SSP limits are noncentral-t and
    genuinely asymmetric: at the 500-year the upper side is about 1.5x the
    lower side. TRANSFORM_SIGMA_ASYMMETRIC extends the same treatment to the
    transform's own scatter -- computed separately from the residuals lying
    above the fitted line versus below it -- rather than assuming by default
    that the cloud of points scatters the same amount on both sides. Each is
    modelled as its own Normal, per EM Sec 4-4.c(1)/4-6.a's convention (whole
    distribution not Normal; each HALF, split at the best estimate, treated as
    one) -- a two-piece normal, still closed form. UNCERTAINTY_REPORT prints
    the combined asymmetry ratio at each AEP; if it gets far from 1 (past
    ASYMMETRY_MC_TRIGGER) the two-piece approximation is being asked to do too
    much and a simple Monte Carlo over the two terms would be the honest next
    step -- EM Sec 4-6 sanctions this too ("if it is not possible to develop
    [an analytical] rating... sensitivity analysis... is an appropriate
    approach").

    PREDICTION vs MEAN LINE. TRANSFORM_UNCERTAINTY_BASIS decides what the
    transform term means:
      "prediction" (default) the full scatter -- "the next flood of this
          magnitude could have any of the shapes we have seen". This is the
          design-relevant question and the wider band.
      "mean"       scatter / sqrt(n) -- "how well is the average relationship
          located". Much narrower, and it answers a question nobody is asking
          of a design curve.

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

# --- synthetic ResSim results ------------------------------------------------
# The routed synthetic members: each was BUILT to an unregulated peak taken off
# the unregulated frequency curve, then routed through ResSim. That makes them
# the only (unreg, reg) pairs that exist above the observed record, which is
# exactly the range the regulated curve has to be drawn through.
SYNTH_RESULTS_CSV = r"../output/synthetic_results.csv"
SHOW_SYNTH_POINTS = True
# Unregulated side: the peak the member was scaled to hit, which is a point on
# the unregulated frequency curve by construction. Regulated side: what ResSim
# routed it to.
# The ROUTED unregulated peak at Castle Rock, from ResSim's Flow-UNREG record.
# This is the right partner for reg_peak: same quantity, same point, so the
# ratio of the two is the reservoir effect and nothing else. The build target
# (target_unreg_peak_cfs) is an unrouted sum of Mossyrock inflow and local, so
# it sits 1-3% above what actually arrives at Castle Rock -- using it would
# credit the reservoir with an attenuation the river reach did.
SYNTH_UNREG_COL = "unreg_peak_routed_cfs"
# Used if the routed record is missing from an older results CSV.
SYNTH_UNREG_COL_FALLBACK = "target_unreg_peak_cfs"
SYNTH_REG_COL = "reg_peak"
# Which source event each member was scaled from. Distinguishing these on the
# plots matters: the spread among synthetics at one magnitude is almost entirely
# a SHAPE effect, not scatter. A sustained event exhausts storage and passes
# through near the 1:1 line, a sharp one is held back hard, and the two land
# 100,000 cfs apart at the same unregulated peak. A single marker colour hides
# that and makes the difference look like noise.
SYNTH_EVENT_COL = "event"
# One marker shape per source event, all in C_SYNTH. Shape rather than colour
# because the colours here are already carrying source (adjusted / WCM_RC /
# synthetic / 2009) and reusing them for event identity would be ambiguous.
# One per source event; the catalog can run to eight or more, and a wrapped
# marker would make two different events look like the same one.
SYNTH_EVENT_MARKERS = ["^", "v", "s", "D", "P", "X", "<", ">", "*", "h"]
# Restrict to one scaling method when both were run.
SYNTH_SCALING_METHOD = "volume_matched"

# --- WHICH POINTS DRAW THE REGULATED LINE ------------------------------------
# Showing a point and fitting a point are different decisions, so they are
# separate settings. Everything switched on above is DRAWN; only what is
# switched on here MOVES THE LINE.
#
#   "adjusted"  the historic adjusted record -- observed peaks corrected by the
#               WCM_RC/Obs_RC difference, screened. Real events, real routing.
#   "wcm_rc"    the raw WCM_RC simulated pairs. These are the PRE-ADJUSTMENT
#               points: no Obs_RC data exists behind them, so no correction was
#               ever applied and their position carries the model's bias. Shown
#               for shape, never used to place the line.
#   "synthetic" the routed synthetic members. The only evidence above the
#               observed record.
FIT_SOURCES = {"adjusted": True, "wcm_rc": False, "synthetic": True}
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
SHOW_POWER_LAW_REFERENCE = False

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
C_SYNTH = "#d68910"        # routed synthetic members
C_2009 = "#117a65"

# The clean figure: regulated and unregulated lines only, no scatter.
MAKE_FINAL_CURVE_PLOT = True

# --- uncertainty on the regulated curve (PART 4 in the docstring) ------------
# The confidence limits HEC-SSP wrote into the frequency table. Read off the
# analysis XML (<ConfidenceLevel>0.05</> and <ConfidenceLevel>0.95</>), i.e. a
# 90% two-sided interval. If the SSP analysis is re-run at a different level
# this MUST be changed to match, or every band below is silently mis-scaled --
# the limits are converted to a sigma by dividing by z at this level.
SSP_CONF_LIMITS = (0.05, 0.95)
FREQ_LOWER_COL = "LowerConf"
FREQ_UPPER_COL = "UpperConf"
# The limits are reported about the COMPUTED curve ("Value"), while the adopted
# unregulated curve is the expected-probability one ("Expected"). The spread is
# therefore taken as a log RATIO to Value and applied to Expected, rather than
# being used as absolute flows.
FREQ_LIMIT_BASE_COL = "Value"
# Cross-check the limits against sqrt(VarianceLog), which is the same standard
# error by a different route. Warn if they disagree by more than this fraction
# -- that means SSP_CONF_LIMITS does not match the analysis that wrote the file.
FREQ_VARIANCE_COL = "VarianceLog"
FREQ_SIGMA_CHECK_TOL = 0.15

# What the FINAL combined band spans. 0.95 (97.5%/2.5%) per the reviewer's
# formula and EM 1110-2-1619 Sec 4-4.c(1): "the range created by the mean plus
# and minus two standard deviations spans 95% of the probability" -- the EM's
# own rounded z=2 convention for a 95% two-sided band; this script uses the
# exact z=1.960 (norm.ppf(0.975)) rather than the rounded 2.
#
# This is the ONLY confidence level that matters for the combination itself.
# Each SOURCE below (frequency, transform) is first reduced to its own proper
# standard deviation, in whatever units/level is natural for that source; the
# two sigmas are combined by root-sum-of-squares (EM eq 4-6); THEN this one z
# is applied once, to the combined sigma. Sources are never combined "at their
# own levels" against each other -- that would only be valid if every source
# happened to share the same z, and 90% (SSP's native reporting) and whatever
# level the transform scatter is estimated at are not guaranteed to match.
UNCERTAINTY_CONF_LEVEL = 0.95
#   "prediction" : full scatter about the transform -- the next flood of this
#                  size could have any observed shape. The design question.
#   "mean"       : scatter / sqrt(n), the uncertainty of the fitted line only.
TRANSFORM_UNCERTAINTY_BASIS = "prediction"
# How the transform's scatter is measured along the curve.
#   "local"    : tricube-weighted residual RMS in the same neighbourhood the
#                LOESS line was drawn from, so the band follows the data. The
#                scatter really does grow with magnitude here -- about 0.05 dex
#                below 80,000 cfs against 0.08 above 150,000 -- because whether
#                a big flood exhausts storage depends on its shape, while a
#                small one is simply held. The default.
#   "constant" : one pooled sd everywhere. Simpler to describe, but it is about
#                25% too wide at the median and 15% too narrow at the top, and
#                it draws the band as a constant ribbon.
TRANSFORM_SIGMA_MODE = "local"
# Neighbourhood for the local scatter. Wider than LOESS_SPAN would over-smooth
# the variance back towards constant; much narrower gets noisy at n=88.
TRANSFORM_SIGMA_SPAN = 0.50
# How the "frequency" (unregulated-curve) term of the combination is
# expressed, before it is RSS'd against the transform term. The reviewer's
# formula is written as (Unreg_bound - Unreg_best) -- a delta in UNREGULATED
# cfs -- added directly to a REGULATED-cfs transform term. That mixes units
# (a regulated flow cannot literally gain "unregulated cfs of uncertainty"),
# and while the two modes track each other reasonably closely over most of
# THIS curve, that is a property of this transform's particular shape (its
# attenuation ratio and log-log slope partly offsetting), not something
# guaranteed for a different basin or a re-fit transform -- see
# freq_term_mode_comparison.txt. "literal" implicitly assumes the transform's
# local slope is 1 (pure pass-through) wherever it skips the conversion,
# which is false through most of this curve (30-37% attenuation). Switched
# to "transform_curve" for that reason -- it is correct by construction
# rather than by coincidence, and costs nothing (the fitted curve is already
# on hand).
#   "literal"         reviewer's formula exactly: freq term = Unreg_bound -
#                      Unreg_best, in raw unregulated cfs. Kept only for
#                      comparison (freq_term_mode_comparison.txt/plots) and
#                      as the literal reading of the reviewer's request.
#   "transform_curve"  push the unregulated bound through the fitted transform
#                      curve first (apply_transform), so the freq term is
#                      already in regulated cfs before combining. Dimensionally
#                      consistent -- ADOPTED. Both modes are always computed
#                      and written to the CSV/report; this only picks which
#                      one is the adopted band vs. the "_alt" comparison
#                      columns.
FREQ_TERM_MODE = "transform_curve"
# Keep the two sides of the frequency interval separate so the noncentral-t
# asymmetry survives the combination. False averages them into one sigma.
COMBINE_ASYMMETRIC = True
# Split the TRANSFORM (LOESS) scatter into an upper-half and lower-half sigma
# too, from the residuals that lie above the fitted line versus below it,
# instead of one pooled number applied to both sides. Tests for a real skew in
# the scatter rather than assuming symmetry the way the old version did.
# False restores one pooled sigma for both sides (still locally-varying by
# magnitude if TRANSFORM_SIGMA_MODE = "local", just not split by sign).
TRANSFORM_SIGMA_ASYMMETRIC = True
# A side needs at least this many EFFECTIVE points (tricube-weighted, not a
# raw count) in the local window before its own sigma is trusted. Below it,
# that side falls back to the pooled sigma at that point, and the fallback is
# counted and reported rather than left silent -- a fallback point looks
# identical to a genuinely symmetric one otherwise.
TRANSFORM_SIDE_MIN_N_EFF = 4.0
# Print the term-by-term table so the band can be audited at each AEP.
UNCERTAINTY_REPORT = True
# Flag if the combined band gets more lopsided than this -- past here the
# two-piece lognormal is being over-asked and Monte Carlo is worth doing.
ASYMMETRY_MC_TRIGGER = 2.0
# Run the Monte Carlo check itself rather than just flagging that it would be
# worth doing. Draws MC_N_DRAWS samples per AEP from the two independent
# sources -- each as a split (two-piece) Normal using the same sigma_lo/hi
# already in unc -- pushes the frequency draws through the actual fitted
# transform curve (not the local-slope approximation), multiplies in the
# transform-side scatter, and compares the empirical 2.5/97.5 percentiles to
# the closed-form band. The closed-form combination (RSS the upper sigmas,
# RSS the lower sigmas) is itself an approximation once two skewed
# distributions are involved; this checks how much that approximation is
# costing, concentrated exactly where ASYMMETRY_MC_TRIGGER says to check it.
# ~50 microseconds per curve evaluation, so 50,000 draws x 16 AEPs runs in
# well under a minute -- turn off only if that is not affordable right now.
RUN_MONTE_CARLO_CHECK = True
MC_N_DRAWS = 50000
MC_SEED = 20260819
# Make the FREQ_TERM_MODE comparison output: two figures (one per mode, each
# alone) plus one figure with both bands overlaid so they can be told apart
# directly, and a plain-text methodology writeup covering both. This is for
# deciding which mode belongs in the memo, not for the memo itself -- nothing
# here touches the adopted band, which is still whichever mode FREQ_TERM_MODE
# is set to.
MAKE_FREQ_MODE_COMPARISON = True

# Hold the regulated UPPER bound at or below the unregulated upper bound at the
# same confidence level. Without this the band goes non-physical at the top: the
# local slope b reaches about 1.5 where the transform bends towards pass-through,
# which amplifies the frequency sigma past the unregulated curve's own sigma, and
# the raw combination ends up claiming a regulated flood larger than the
# unregulated flood it was routed from. The central line is already clipped at
# 1:1 by CLIP_TO_UNREG; this applies the same physics to the band.
#
# Only above BAND_CLIP_MIN_CFS. Below it reg > unreg is real and documented --
# minimum release and refill drawdown put more water in the river than nature
# would -- so clipping there would be wrong. Same threshold and same reasoning
# as the reg-over-unreg screen in #Adjusted_Peak_Record.py.
CLIP_BAND_TO_UNREG = True
BAND_CLIP_MIN_CFS = 60000.0

# The adopted report figure. The 2009 curve is a diagnostic comparison and is
# deliberately left off it.
MAKE_FINAL_UNCERTAINTY_PLOT = True
FINAL_SHOW_2009 = False
# Draw the unregulated confidence band as well, for context.
FINAL_SHOW_UNREG_BAND = True
# Title. One line, no explanation underneath -- this figure goes in a memo
# where the method is in the text beside it.
FINAL_TITLE = "Castle Rock peak flow frequency, regulated and unregulated"
# The dotted line and label marking where the transform stops being supported
# by data. Off for the report figure; the caveat belongs in the memo text,
# where it can be stated properly rather than in eight-point grey. The other
# plots still carry it -- see plot_frequency.
FINAL_SHOW_SUPPORT_MARKER = False
# The box in the lower right restating the band formula. Off for the same
# reason. Switch it on when the figure has to travel on its own.
FINAL_SHOW_FORMULA_NOTE = False
# Frequency axis for THIS figure only. The global AEP_LIMITS starts at 0.999,
# which leaves a wide empty strip on the left because no curve is plotted
# beyond 0.99. Kept separate so the diagnostic plots are not moved with it.
FINAL_AEP_LIMITS = (0.99, 0.0005)

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
        lx, ly = np.log10(x), np.log10(y)
        fit["lx"], fit["ly"] = lx, ly
        fit["resid"] = ly - np.log10(apply_power_law(power, x))
        fit["se_dex"] = power["se_dex"]
        fit["r2"] = power["r2"]
        return fit

    lx, ly = np.log10(x), np.log10(y)
    fit["lx"], fit["ly"] = lx, ly
    fit["span"] = LOESS_SPAN
    predicted = np.array([loess_at(lx, ly, v, LOESS_SPAN) for v in lx])
    resid = ly - predicted
    fit["resid"] = resid
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


def load_synthetic_points(csv_path):
    """Routed synthetic members as (unreg, reg) pairs.

    The unregulated side is the peak the member was BUILT to hit, not a peak
    read back out of a simulation, so it is exact by construction and sits on
    the unregulated frequency curve. The regulated side is what ResSim did with
    it.
    """
    if not (SHOW_SYNTH_POINTS and os.path.isfile(csv_path)):
        if SHOW_SYNTH_POINTS:
            print("Synthetic : not found: %s -- points not drawn" % csv_path)
        return None
    table = pd.read_csv(csv_path)
    if SYNTH_SCALING_METHOD and "scaling_method" in table.columns:
        table = table[table["scaling_method"] == SYNTH_SCALING_METHOD]
    global SYNTH_UNREG_COL
    if SYNTH_UNREG_COL not in table.columns or table[SYNTH_UNREG_COL].isna().all():
        print("Synthetic : no '%s' -- falling back to '%s'. That column is the"
              % (SYNTH_UNREG_COL, SYNTH_UNREG_COL_FALLBACK))
        print("            UNROUTED build target, so the pairs will overstate")
        print("            attenuation by the routing loss. Re-run")
        print("            #Synthetic_Diagnostics.py to get the routed peak.")
        SYNTH_UNREG_COL = SYNTH_UNREG_COL_FALLBACK
    for col in (SYNTH_UNREG_COL, SYNTH_REG_COL):
        if col not in table.columns:
            print("Synthetic : %s has no '%s' column -- points not drawn"
                  % (os.path.basename(csv_path), col))
            return None
    table = table.dropna(subset=[SYNTH_UNREG_COL, SYNTH_REG_COL]).copy()
    table = table[(table[SYNTH_UNREG_COL] > 0) & (table[SYNTH_REG_COL] > 0)]
    return table.reset_index(drop=True)


def assemble_fit_points(data, wcm, synth):
    """Stack the point sets that are switched on in FIT_SOURCES.

    Returns the arrays and a per-source count, so the log states exactly what
    moved the line rather than leaving it to be inferred from the plot.
    """
    xs, ys, used = [], [], {}
    if FIT_SOURCES.get("adjusted", False) and data is not None and len(data):
        xs.append(np.asarray(data[UNREG_COL].values, dtype=float))
        ys.append(np.asarray(data["reg_peak"].values, dtype=float))
        used["adjusted"] = len(data)
    if FIT_SOURCES.get("wcm_rc", False) and wcm is not None and len(wcm):
        xs.append(np.asarray(wcm["unreg_peak"].values, dtype=float))
        ys.append(np.asarray(wcm["reg_peak"].values, dtype=float))
        used["wcm_rc"] = len(wcm)
    if FIT_SOURCES.get("synthetic", False) and synth is not None and len(synth):
        xs.append(np.asarray(synth[SYNTH_UNREG_COL].values, dtype=float))
        ys.append(np.asarray(synth[SYNTH_REG_COL].values, dtype=float))
        used["synthetic"] = len(synth)
    if not xs:
        raise SystemExit("FIT_SOURCES leaves nothing to fit -- switch at least "
                         "one point set on.")
    return np.concatenate(xs), np.concatenate(ys), used


def synth_by_event(synth):
    """Yield (event, marker, rows) so each source event is its own plot series.

    Events keep the order they appear in the results CSV, so the marker assigned
    to an event is stable across the scatter, the frequency plot and the log.
    """
    if SYNTH_EVENT_COL not in synth.columns:
        yield "synthetic", SYNTH_EVENT_MARKERS[0], synth
        return
    for i, event in enumerate(dict.fromkeys(synth[SYNTH_EVENT_COL])):
        marker = SYNTH_EVENT_MARKERS[i % len(SYNTH_EVENT_MARKERS)]
        yield str(event), marker, synth[synth[SYNTH_EVENT_COL] == event]


def report_synth_by_event(synth):
    """What each source event does to a flood, in numbers rather than by eye.

    The attenuation ratio is the whole story: near 1.0 means the project passed
    the flood through, low means it held it back. Printing it per event says
    which markers are which without anyone having to read the plot.
    """
    if synth is None or not len(synth) or SYNTH_EVENT_COL not in synth.columns:
        return
    print("\nSYNTHETICS BY SOURCE EVENT  (reg / unreg, 1.00 = full pass-through)")
    for event, marker, group in synth_by_event(synth):
        ratio = group[SYNTH_REG_COL] / group[SYNTH_UNREG_COL]
        print("   '%s'  %-9s n=%d   attenuation %.2f to %.2f   reg %s to %s cfs"
              % (marker, event, len(group), ratio.min(), ratio.max(),
                 format(int(group[SYNTH_REG_COL].min()), ","),
                 format(int(group[SYNTH_REG_COL].max()), ",")))
    # Median, not max: the question is which EVENT passes floods through, not
    # which single member happened to be largest. One member of a well-attenuated
    # event can reach 1:1 at the top magnitude without that event being the
    # pass-through case.
    top = max(synth_by_event(synth),
              key=lambda g: (g[2][SYNTH_REG_COL] / g[2][SYNTH_UNREG_COL]).median())
    print("   Sitting on the 1:1 line: %s -- storage is exhausted at every"
          % top[0])
    print("   magnitude and the flood passes through, so it sets the upper end")
    print("   of the regulated curve.")


def annotate_points(ax, x, y, labels):
    """Callouts placed alternately above and below to reduce collisions."""
    for i, (xi, yi, text) in enumerate(zip(x, y, labels)):
        dy = 14 if i % 2 == 0 else -18
        ax.annotate(text, (xi, yi), xytext=(10, dy), textcoords="offset points",
                    fontsize=8, color=C_REG,
                    arrowprops=dict(arrowstyle="-", color=C_REG, lw=0.7))


def plot_scatter(data, fit, wcm, synth, stem):
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
                       label="Unadjusted%s)"
                             % ("" if FIT_SOURCES.get("wcm_rc") else ", (not fitted"))
            lo += [wcm["unreg_peak"].min(), wcm["reg_peak"].min()]
            hi += [wcm["unreg_peak"].max(), wcm["reg_peak"].max()]
        if synth is not None and len(synth):
            for event, marker, group in synth_by_event(synth):
                ax.scatter(group[SYNTH_UNREG_COL], group[SYNTH_REG_COL], s=70,
                           marker=marker, facecolor=C_SYNTH, edgecolor="0.2",
                           lw=0.6, zorder=4,
                           label="Synthetic from %s (n=%d)" % (event, len(group)))
            lo += [synth[SYNTH_UNREG_COL].min(), synth[SYNTH_REG_COL].min()]
            hi += [synth[SYNTH_UNREG_COL].max(), synth[SYNTH_REG_COL].max()]

        ax.scatter(x, y, s=46, facecolor=C_UNREG, edgecolor="0.25", lw=0.6,
                   zorder=3, label="Adjusted")
        lim = [min(lo) * 0.85, max(hi) * 1.15]
        ax.plot(lim, lim, color="k", lw=1.2, ls="--", zorder=2,
                label="1:1")

        xs = np.geomspace(max(lim[0], 1.0), lim[1], 300)
        centre = apply_transform(fit, xs)
        ax.plot(xs, centre, color=C_REG, lw=1.8, zorder=4,
                label=transform_label(fit))
        # band = 10 ** transform_sigma_dex(fit, xs)
        # ax.fill_between(xs, centre / band, centre * band, color=C_REG,
        #                 alpha=0.12, zorder=1,
        #                 label="+/- 1 sigma of the scatter (%s)"
        #                       % TRANSFORM_SIGMA_MODE)
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


def plot_frequency(freq, data, fit, wcm, synth, table_2009, stem):
    """Unregulated and inferred regulated frequency curves, SSP idiom."""
    fig, ax = plt.subplots(figsize=(11.5, 8.8))

    aep_values = freq["AEP"].values
    z = stats.norm.ppf(1.0 - aep_values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    reg_curve = apply_transform(fit, unreg_curve)
    band = 10 ** transform_sigma_dex(fit, unreg_curve)

    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.2, zorder=4,
            label="Unregulated peak (%s curve)" % FREQ_VALUE_COL.lower())
    ax.plot(z, reg_curve, color=C_REG, lw=2.2, zorder=4,
            label="Regulated peak, inferred (%s)"
                  % ("LOESS" if fit["method"] == "loess"
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
        # --- routed synthetics, at the AEP their unreg target came from ------
        if synth is not None and len(synth):
            aep_s = aep_from_unreg_curve(synth[SYNTH_UNREG_COL].values,
                                         unreg_curve, aep_values)
            zs = stats.norm.ppf(1.0 - np.clip(aep_s, 1e-6, 1 - 1e-6))
            synth = synth.copy()
            synth["_z"] = zs
            for event, marker, group in synth_by_event(synth):
                ax.plot(group["_z"], group[SYNTH_REG_COL].values, ls="none",
                        marker=marker, ms=8, mfc=C_SYNTH, mec="0.2", mew=0.6,
                        zorder=7, label="Synthetic regulated, from %s" % event)
            ax.plot(zs, synth[SYNTH_UNREG_COL].values, ls="none", marker="_",
                    ms=9, mfc="none", mec=C_SYNTH, mew=1.2, zorder=6,
                    label="Synthetic unregulated (routed)")
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
    # supported = unreg_curve <= fit["x_max"]
    # if supported.any() and not supported.all():
    #     z_edge = z[supported].max()
    #     ax.axvline(z_edge, color="0.35", lw=1.0, ls=":", zorder=2)
    #     # ax.text(z_edge, FLOW_LIMITS[1] * 0.92,
    #     #         "  transform supported to here\n  (unreg %s cfs)"
    #     #         % format(int(fit["x_max"]), ","),
    #     #         fontsize=8, color="0.3", va="top")

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
                 % ("(LOESS)" if fit["method"] == "loess"
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


def frequency_sigma_dex(freq):
    """Log10 standard error of the unregulated quantile, per side, per AEP.

    The SSP confidence limits are noncentral-t and asymmetric, so each side is
    converted to its own sigma by dividing by z at the level SSP used. Taken as
    a RATIO to FREQ_LIMIT_BASE_COL (the computed curve the limits belong to) so
    it can be applied to the adopted Expected curve without mixing the two.

    Cross-checked against sqrt(VarianceLog), which is the same standard error
    by a different route: if the mean of the two sides disagrees with it, the
    declared SSP_CONF_LIMITS does not match the analysis that wrote the table
    and every band downstream is mis-scaled.
    """
    z = stats.norm.ppf(SSP_CONF_LIMITS[1])
    base = freq[FREQ_LIMIT_BASE_COL].values.astype(float)
    upper = freq[FREQ_UPPER_COL].values.astype(float)
    lower = freq[FREQ_LOWER_COL].values.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_hi = np.log10(upper / base) / z
        sigma_lo = np.log10(base / lower) / z
    sigma_hi = np.where(np.isfinite(sigma_hi) & (sigma_hi > 0), sigma_hi, np.nan)
    sigma_lo = np.where(np.isfinite(sigma_lo) & (sigma_lo > 0), sigma_lo, np.nan)

    note = ""
    if FREQ_VARIANCE_COL in freq.columns:
        declared = np.sqrt(freq[FREQ_VARIANCE_COL].values.astype(float))
        mean_side = 0.5 * (sigma_hi + sigma_lo)
        good = np.isfinite(declared) & (declared > 0) & np.isfinite(mean_side)
        if good.any():
            rel = np.abs(mean_side[good] / declared[good] - 1.0)
            note = ("limits vs sqrt(%s): %.1f%% to %.1f%% apart"
                    % (FREQ_VARIANCE_COL, 100 * rel.min(), 100 * rel.max()))
            if rel.max() > FREQ_SIGMA_CHECK_TOL:
                print("WARNING the SSP confidence limits do not match "
                      "sqrt(%s) to within %.0f%%." % (FREQ_VARIANCE_COL,
                                                      100 * FREQ_SIGMA_CHECK_TOL))
                print("        %s" % note)
                print("        SSP_CONF_LIMITS is set to %s -- check it against "
                      "the" % (SSP_CONF_LIMITS,))
                print("        ConfidenceLevel entries in the SSP analysis XML.")
    return sigma_hi, sigma_lo, note


def transform_log_slope(fit, x_eval):
    """Local slope b = dlog10(reg)/dlog10(unreg) of the transform.

    Taken by finite difference off a dense evaluation of the SAME apply_transform
    the curve is drawn with, so the monotonic enforcement and the 1:1 clip are
    both inside the slope rather than being ignored by it. Evaluating the LOESS
    directly at two nearby points would miss them.
    """
    x_eval = np.asarray(x_eval, dtype=float)
    grid = np.geomspace(max(x_eval.min() * 0.75, 1.0), x_eval.max() * 1.3, 400)
    y = apply_transform(fit, grid)
    lg, ly = np.log10(grid), np.log10(np.clip(y, 1e-9, None))
    slope = np.gradient(ly, lg)
    return np.interp(np.log10(x_eval), lg, slope)


def local_sigma_at(lx, resid, x0, span):
    """Tricube-weighted RMS of the residuals near x0, in log10 space.

    The variance companion to loess_at: the same neighbourhood and the same
    weights, so the scatter is measured where the line was drawn. Corrected by
    the effective sample size of the weights rather than a raw count, since the
    tricube kernel means the far points in the window barely contribute.
    """
    n = len(lx)
    k = int(np.ceil(span * n))
    k = max(8, min(k, n))
    dist = np.abs(lx - x0)
    near = np.argsort(dist)[:k]
    d = dist[near]
    d_max = d.max()
    weights = np.ones(k) if d_max <= 0 else (1.0 - (d / d_max) ** 3) ** 3
    weights = np.clip(weights, 1e-8, None)
    r = resid[near]
    w_sum = weights.sum()
    n_eff = w_sum ** 2 / np.sum(weights ** 2)
    var = np.sum(weights * r * r) / w_sum
    return float(np.sqrt(var * n_eff / max(n_eff - 2.0, 1.0)))


def transform_sigma_dex(fit, x_eval):
    """The transform's log10 sigma at each unregulated flow.

    Returns an array, because the scatter is NOT constant along the curve. The
    reservoir is predictable at small floods -- it simply holds them -- and much
    less so at large ones, where the answer depends on whether that particular
    flood exhausts storage. Measured on this dataset the residual sd runs about
    0.05 dex below 80,000 cfs and 0.08 dex above 150,000, a 1.5x spread.

    A single pooled se_dex splits the difference and is therefore wrong at both
    ends -- roughly 25% too wide at the median and 15% too narrow at the top,
    which is exactly where the curve is being used. It also makes the band look
    like a constant ribbon, because it is one.

    TRANSFORM_SIGMA_MODE = "constant" restores the pooled value.
    """
    x_eval = np.atleast_1d(np.asarray(x_eval, dtype=float))
    if TRANSFORM_SIGMA_MODE == "constant" or "resid" not in fit:
        sigma = np.full(len(x_eval), fit["se_dex"], dtype=float)
    else:
        lx, resid = fit["lx"], fit["resid"]
        sigma = np.array([local_sigma_at(lx, resid, v, TRANSFORM_SIGMA_SPAN)
                          for v in np.log10(np.clip(x_eval, 1e-6, None))])
    if TRANSFORM_UNCERTAINTY_BASIS == "mean":
        sigma = sigma / np.sqrt(max(fit["n"], 1))
    return sigma


def local_sigma_at_side(lx, resid, x0, span, side, min_n_eff):
    """local_sigma_at, restricted to residuals on one side of the fitted line.

    side="hi" uses only resid > 0 (points that lie above the line), "lo" only
    resid < 0. Same neighbourhood and same tricube weights as local_sigma_at --
    only the residuals actually summed differ -- so a side sigma is comparable
    to the pooled one, not a different statistic in disguise.

    Returns None if that side does not have enough EFFECTIVE weight nearby
    (min_n_eff) to trust its own number; the caller falls back to the pooled
    sigma at that point and counts the fallback.
    """
    n = len(lx)
    k = int(np.ceil(span * n))
    k = max(8, min(k, n))
    dist = np.abs(lx - x0)
    near = np.argsort(dist)[:k]
    d = dist[near]
    d_max = d.max()
    weights = np.ones(k) if d_max <= 0 else (1.0 - (d / d_max) ** 3) ** 3
    weights = np.clip(weights, 1e-8, None)
    r = resid[near]
    mask = (r > 0) if side == "hi" else (r < 0)
    if not mask.any():
        return None
    w_sub, r_sub = weights[mask], r[mask]
    w_sum = w_sub.sum()
    n_eff = w_sum ** 2 / np.sum(w_sub ** 2)
    if n_eff < min_n_eff:
        return None
    var = np.sum(w_sub * r_sub * r_sub) / w_sum
    return float(np.sqrt(var * n_eff / max(n_eff - 2.0, 1.0)))


def transform_sigma_dex_split(fit, x_eval):
    """Upper- and lower-half log10 sigma of the transform's scatter.

    Same idea as transform_sigma_dex, but the residuals lying ABOVE the fitted
    line and those lying BELOW it are summarised separately, so a genuine skew
    in the scatter (more room on one side of the line than the other) survives
    into the band instead of being pooled into one symmetric number.

    TRANSFORM_SIGMA_ASYMMETRIC = False returns the pooled sigma on both sides,
    i.e. the old symmetric behaviour. Returns (sigma_hi, sigma_lo, n_fallback,
    n_sides) -- the last two so the caller can report how often a side had to
    fall back to the pooled value for lack of nearby data on that side alone.
    """
    x_eval = np.atleast_1d(np.asarray(x_eval, dtype=float))
    pooled = transform_sigma_dex(fit, x_eval)
    if (not TRANSFORM_SIGMA_ASYMMETRIC or TRANSFORM_SIGMA_MODE == "constant"
            or "resid" not in fit):
        return pooled.copy(), pooled.copy(), 0, 0

    lx, resid = fit["lx"], fit["resid"]
    log_x = np.log10(np.clip(x_eval, 1e-6, None))
    sigma_hi = np.empty(len(x_eval))
    sigma_lo = np.empty(len(x_eval))
    n_fallback = 0
    for i, v in enumerate(log_x):
        hi = local_sigma_at_side(lx, resid, v, TRANSFORM_SIGMA_SPAN, "hi",
                                 TRANSFORM_SIDE_MIN_N_EFF)
        lo = local_sigma_at_side(lx, resid, v, TRANSFORM_SIGMA_SPAN, "lo",
                                 TRANSFORM_SIDE_MIN_N_EFF)
        if hi is None:
            hi, n_fallback = pooled[i], n_fallback + 1
        if lo is None:
            lo, n_fallback = pooled[i], n_fallback + 1
        sigma_hi[i], sigma_lo[i] = hi, lo
    if TRANSFORM_UNCERTAINTY_BASIS == "mean":
        sigma_hi = sigma_hi / np.sqrt(max(fit["n"], 1))
        sigma_lo = sigma_lo / np.sqrt(max(fit["n"], 1))
    return sigma_hi, sigma_lo, n_fallback, 2 * len(x_eval)


def combine_uncertainty(freq, fit, reg_curve):
    """Regulated curve confidence band, matching the reviewer's formula:

        Upper = RegBest + sqrt((Unreg_hi - Unreg_best)^2 + (Transform_hi - Transform_best)^2)
        Lower = RegBest - sqrt((Unreg_best - Unreg_lo)^2 + (Transform_best - Transform_lo)^2)

    at UNCERTAINTY_CONF_LEVEL (0.95 two-sided -> 97.5%/2.5%, per the reviewer's
    own example), with each side's sigma kept separate (EM 1110-2-1619 sec
    4-4.c(1)/4-6a: upper half and lower half both treated as Normal, each with
    its own standard deviation, rather than assuming one symmetric sigma).

    The frequency term is computed BOTH ways every run (FREQ_TERM_MODE picks
    which one is used for the adopted band; both are returned for
    comparison/reporting):
      "literal"          Unreg_hi - Unreg_best, in raw unregulated cfs -- the
                          reviewer's formula exactly.
      "transform_curve"  the unregulated bound pushed through the same fitted
                          transform curve as the central estimate, so the term
                          is in regulated cfs before combining -- dimensionally
                          consistent, kept for comparison only unless selected.
    The transform term is Transform_hi - Transform_best, i.e. the LOESS fit's
    own confidence band at the unregulated BEST ESTIMATE, already in regulated
    cfs either way.

    Returns the band and every term that went into it, so the CSV can show the
    working and the report can print it AEP by AEP.
    """
    sigma_freq_hi, sigma_freq_lo, note = frequency_sigma_dex(freq)
    if not COMBINE_ASYMMETRIC:
        mean_side = np.nanmean(np.vstack([sigma_freq_hi, sigma_freq_lo]), axis=0)
        sigma_freq_hi = sigma_freq_lo = mean_side

    unreg_best = freq[FREQ_VALUE_COL].values.astype(float)
    z = stats.norm.ppf(0.5 + UNCERTAINTY_CONF_LEVEL / 2.0)

    # --- frequency term, both ways --------------------------------------
    unreg_upper = unreg_best * 10.0 ** (z * sigma_freq_hi)
    unreg_lower = unreg_best * 10.0 ** (-z * sigma_freq_lo)

    # "literal": the reviewer's formula exactly -- delta in unregulated cfs.
    freq_term_hi_literal = np.maximum(unreg_upper - unreg_best, 0.0)
    freq_term_lo_literal = np.maximum(unreg_best - unreg_lower, 0.0)

    # "transform_curve": same delta, but pushed through the fitted curve
    # first so it is in regulated cfs before combining. apply_transform is
    # monotonic non-decreasing (ENFORCE_MONOTONIC), so these are >= 0
    # already; the floor is only a guard against float noise on a flat
    # stretch of the curve, not a correction to anything real.
    reg_from_unreg_upper = apply_transform(fit, unreg_upper)
    reg_from_unreg_lower = apply_transform(fit, unreg_lower)
    freq_term_hi_curve = np.maximum(reg_from_unreg_upper - reg_curve, 0.0)
    freq_term_lo_curve = np.maximum(reg_curve - reg_from_unreg_lower, 0.0)

    if FREQ_TERM_MODE == "transform_curve":
        freq_term_hi, freq_term_lo = freq_term_hi_curve, freq_term_lo_curve
    else:
        freq_term_hi, freq_term_lo = freq_term_hi_literal, freq_term_lo_literal

    # --- transform term, held at the unregulated BEST ESTIMATE --------------
    sigma_t_hi, sigma_t_lo, n_fallback, n_sides = transform_sigma_dex_split(
        fit, unreg_best)
    transform_term_hi = reg_curve * (10.0 ** (z * sigma_t_hi) - 1.0)
    transform_term_lo = reg_curve * (1.0 - 10.0 ** (-z * sigma_t_lo))

    # --- combine per side, additively, in cfs -- eq 4-6 ----------------------
    delta_hi = np.sqrt(freq_term_hi ** 2 + transform_term_hi ** 2)
    delta_lo = np.sqrt(freq_term_lo ** 2 + transform_term_lo ** 2)
    reg_upper = reg_curve + delta_hi
    reg_lower = np.maximum(reg_curve - delta_lo, 0.0)

    # the other mode's band too, purely for comparison/reporting
    delta_hi_alt = np.sqrt(freq_term_hi_curve ** 2 + transform_term_hi ** 2) \
        if FREQ_TERM_MODE != "transform_curve" \
        else np.sqrt(freq_term_hi_literal ** 2 + transform_term_hi ** 2)
    delta_lo_alt = np.sqrt(freq_term_lo_curve ** 2 + transform_term_lo ** 2) \
        if FREQ_TERM_MODE != "transform_curve" \
        else np.sqrt(freq_term_lo_literal ** 2 + transform_term_lo ** 2)

    # A regulated flood cannot be bigger than the unregulated flood it came
    # from -- see CLIP_BAND_TO_UNREG. Applied above BAND_CLIP_MIN_CFS only.
    clipped = np.zeros(len(reg_upper), dtype=bool)
    if CLIP_BAND_TO_UNREG:
        bites = (reg_upper > unreg_upper) & (unreg_best >= BAND_CLIP_MIN_CFS)
        reg_upper = np.where(bites, unreg_upper, reg_upper)
        clipped = bites

    # slope is no longer used in the calculation -- kept only to explain, in
    # the report, why the frequency term grows disproportionately at the top.
    slope = transform_log_slope(fit, unreg_best)
    with np.errstate(divide="ignore", invalid="ignore"):
        asymmetry = np.where(delta_lo > 0, delta_hi / delta_lo, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        asymmetry_alt = np.where(delta_lo_alt > 0, delta_hi_alt / delta_lo_alt, np.nan)
    alt_mode = "transform_curve" if FREQ_TERM_MODE != "transform_curve" else "literal"
    reg_upper_alt = reg_curve + delta_hi_alt
    reg_lower_alt = np.maximum(reg_curve - delta_lo_alt, 0.0)
    if CLIP_BAND_TO_UNREG:
        bites_alt = (reg_upper_alt > unreg_upper) & (unreg_best >= BAND_CLIP_MIN_CFS)
        reg_upper_alt = np.where(bites_alt, unreg_upper, reg_upper_alt)

    return {
        "reg_upper_clipped": clipped,
        "freq_term_mode": FREQ_TERM_MODE,
        "sigma_freq_hi": sigma_freq_hi, "sigma_freq_lo": sigma_freq_lo,
        "sigma_transform_hi": sigma_t_hi, "sigma_transform_lo": sigma_t_lo,
        "transform_sigma_fallback": n_fallback, "transform_sigma_sides": n_sides,
        "freq_term_hi_cfs": freq_term_hi, "freq_term_lo_cfs": freq_term_lo,
        "transform_term_hi_cfs": transform_term_hi,
        "transform_term_lo_cfs": transform_term_lo,
        "delta_hi_cfs": delta_hi, "delta_lo_cfs": delta_lo,
        "slope": slope, "z": z, "note": note,
        "reg_upper": reg_upper, "reg_lower": reg_lower,
        # the unregulated band at the same level, for context on the plot
        "unreg_upper": unreg_upper, "unreg_lower": unreg_lower,
        "asymmetry": asymmetry,
        # --- the OTHER FREQ_TERM_MODE, computed alongside for comparison ---
        "alt_freq_term_mode": alt_mode,
        "freq_term_hi_literal_cfs": freq_term_hi_literal,
        "freq_term_lo_literal_cfs": freq_term_lo_literal,
        "freq_term_hi_curve_cfs": freq_term_hi_curve,
        "freq_term_lo_curve_cfs": freq_term_lo_curve,
        "reg_upper_alt": reg_upper_alt, "reg_lower_alt": reg_lower_alt,
        "delta_hi_alt_cfs": delta_hi_alt, "delta_lo_alt_cfs": delta_lo_alt,
        "asymmetry_alt": asymmetry_alt,
    }


def report_uncertainty(freq, unc, fit):
    """Print the band term by term so it can be audited rather than trusted."""
    if not UNCERTAINTY_REPORT:
        return
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))
    z = unc["z"]
    print("\nREGULATED CURVE UNCERTAINTY  (%d%% band -> z=%.3f, reviewer's formula)"
          % (pct, z))
    print("   Upper = RegBest + sqrt(freq_term_hi^2 + transform_term_hi^2)")
    print("   Lower = RegBest - sqrt(freq_term_lo^2 + transform_term_lo^2)")
    print("   FREQ_TERM_MODE = %r (adopted band below); alt = %r (comparison "
          "columns *_alt)" % (unc["freq_term_mode"], unc["alt_freq_term_mode"]))
    st_hi = np.atleast_1d(unc["sigma_transform_hi"])
    st_lo = np.atleast_1d(unc["sigma_transform_lo"])
    print("   transform sigma: hi %.4f-%.4f dex, lo %.4f-%.4f dex (%s, %s%s)"
          % (np.nanmin(st_hi), np.nanmax(st_hi), np.nanmin(st_lo), np.nanmax(st_lo),
             TRANSFORM_SIGMA_MODE, TRANSFORM_UNCERTAINTY_BASIS,
             ", scatter/sqrt(%d)" % fit["n"]
             if TRANSFORM_UNCERTAINTY_BASIS == "mean" else ", full scatter"))
    if unc["transform_sigma_fallback"]:
        print("   %d of %d side-evaluations fell back to the pooled transform "
              "sigma (too few effective points on that side)"
              % (unc["transform_sigma_fallback"], unc["transform_sigma_sides"]))
    print("   frequency term : from the SSP %g/%g limits%s"
          % (SSP_CONF_LIMITS[0], SSP_CONF_LIMITS[1],
             " -- %s" % unc["note"] if unc["note"] else ""))
    show = pd.DataFrame({
        "AEP": freq["AEP"].values,
        "sig_freq_lo": unc["sigma_freq_lo"], "sig_freq_hi": unc["sigma_freq_hi"],
        "freq_lo_cfs": unc["freq_term_lo_cfs"], "freq_hi_cfs": unc["freq_term_hi_cfs"],
        "trans_lo_cfs": unc["transform_term_lo_cfs"],
        "trans_hi_cfs": unc["transform_term_hi_cfs"],
        "reg_lower": unc["reg_lower"], "reg_upper": unc["reg_upper"],
        "reg_upper_alt": unc["reg_upper_alt"],
        "asym": unc["asymmetry"],
    })
    show = show[show["AEP"].isin([0.5, 0.1, 0.02, 0.01, 0.005, 0.002, 0.001])]
    print(show.round({"AEP": 4, "sig_freq_lo": 4, "sig_freq_hi": 4,
                      "freq_lo_cfs": 0, "freq_hi_cfs": 0, "trans_lo_cfs": 0,
                      "trans_hi_cfs": 0, "reg_lower": 0, "reg_upper": 0,
                      "reg_upper_alt": 0, "asym": 3}).to_string(index=False))
    n_clip = int(unc["reg_upper_clipped"].sum())
    if n_clip:
        aeps = freq["AEP"].values[unc["reg_upper_clipped"]]
        print("   upper bound held at the unregulated upper bound for %d AEP(s): %s"
              % (n_clip, ", ".join("%g" % a for a in aeps)))
        print("   (b reaches %.2f there, which would otherwise push the regulated"
              % np.nanmax(unc["slope"]))
        print("   band above the unregulated flood it was routed from.)")
    worst = np.nanmax(np.abs(np.log(unc["asymmetry"])))
    worst_ratio = float(np.exp(worst))
    print("   most lopsided AEP: %.2fx between the two sides" % worst_ratio)
    if worst_ratio > ASYMMETRY_MC_TRIGGER:
        print("   *** past ASYMMETRY_MC_TRIGGER (%.1fx). The two-piece lognormal"
              % ASYMMETRY_MC_TRIGGER)
        print("   is being over-asked%s" %
              (" -- see monte_carlo_check.txt" if RUN_MONTE_CARLO_CHECK
               else " -- run monte_carlo_check() instead of trusting this band"))
    else:
        print("   within ASYMMETRY_MC_TRIGGER (%.1fx), so the closed-form"
              % ASYMMETRY_MC_TRIGGER)
        print("   two-piece combination is adequate -- no Monte Carlo needed.")


def sample_split_normal(sigma_lo, sigma_hi, n, rng):
    """n draws from a two-piece (split) Normal: Half-Normal(sigma_lo) below
    zero, Half-Normal(sigma_hi) above, continuous at the join.

    The mixture weight on each side has to be proportional to that side's
    sigma -- p_lo = sigma_lo / (sigma_lo + sigma_hi) -- for the density to
    actually be continuous at zero. A 50/50 split would put a visible step
    in the density at the best estimate whenever sigma_lo != sigma_hi.
    """
    sigma_lo = max(float(sigma_lo), 1e-12)
    sigma_hi = max(float(sigma_hi), 1e-12)
    p_lo = sigma_lo / (sigma_lo + sigma_hi)
    is_lo = rng.random(n) < p_lo
    out = np.empty(n)
    out[is_lo] = -np.abs(rng.normal(0.0, sigma_lo, int(is_lo.sum())))
    out[~is_lo] = np.abs(rng.normal(0.0, sigma_hi, int((~is_lo).sum())))
    return out


def monte_carlo_check(freq, fit, unc, reg_curve, n_draws=MC_N_DRAWS, seed=MC_SEED):
    """Simulate the combined band instead of RSS-ing the two sigmas, and
    compare. See RUN_MONTE_CARLO_CHECK for what this is checking and why.

    Per AEP: draw n_draws frequency-side offsets from a split-Normal using
    sigma_freq_lo/hi, push them through the SAME fitted transform curve used
    everywhere else (apply_transform -- monotonic and 1:1-clipped exactly
    like the deterministic curve), then draw transform-side offsets from a
    split-Normal using sigma_transform_lo/hi (evaluated at Unreg_best, the
    same fixed point the closed form uses -- this isolates the question of
    whether RSS-combining the two sources is accurate, not whether transform
    scatter should vary with the sampled draw, which the closed form does
    not attempt either).

    The reg <= unreg physical floor is applied per draw, but ONLY above
    BAND_CLIP_MIN_CFS -- matching CLIP_BAND_TO_UNREG's own condition, not
    applying it unconditionally. Below that threshold reg > unreg is real
    (minimum release and refill drawdown put more water in the river than
    nature would), so clipping every draw there would manufacture a fake
    ceiling on the low-flow upper tail and understate that part of the band
    for a reason that has nothing to do with the two-piece-lognormal
    question this check exists to answer.

    Returns a DataFrame with the closed-form and MC bounds side by side.
    """
    rng = np.random.default_rng(seed)
    unreg_best = freq[FREQ_VALUE_COL].values.astype(float)
    aep = freq["AEP"].values
    n = len(aep)

    mc_lower = np.empty(n)
    mc_upper = np.empty(n)
    mc_median = np.empty(n)

    for i in range(n):
        freq_dex = sample_split_normal(unc["sigma_freq_lo"][i],
                                       unc["sigma_freq_hi"][i], n_draws, rng)
        unreg_draws = unreg_best[i] * 10.0 ** freq_dex
        reg_from_curve = apply_transform(fit, unreg_draws)

        trans_dex = sample_split_normal(unc["sigma_transform_lo"][i],
                                        unc["sigma_transform_hi"][i],
                                        n_draws, rng)
        reg_draws = reg_from_curve * 10.0 ** trans_dex
        if CLIP_BAND_TO_UNREG and unreg_best[i] >= BAND_CLIP_MIN_CFS:
            reg_draws = np.minimum(reg_draws, unreg_draws)

        mc_lower[i] = np.percentile(reg_draws, 2.5)
        mc_upper[i] = np.percentile(reg_draws, 97.5)
        mc_median[i] = np.percentile(reg_draws, 50.0)

    out = pd.DataFrame({
        "AEP": aep,
        "reg_best_cfs": reg_curve,
        "closed_lower_cfs": unc["reg_lower"],
        "closed_upper_cfs": unc["reg_upper"],
        "mc_lower_cfs": mc_lower,
        "mc_median_cfs": mc_median,
        "mc_upper_cfs": mc_upper,
    })
    with np.errstate(divide="ignore", invalid="ignore"):
        out["lower_diff_pct"] = 100.0 * (out["closed_lower_cfs"] -
                                         out["mc_lower_cfs"]) / out["mc_lower_cfs"]
        out["upper_diff_pct"] = 100.0 * (out["closed_upper_cfs"] -
                                         out["mc_upper_cfs"]) / out["mc_upper_cfs"]
    return out


def plot_monte_carlo_check(freq, fit, reg_curve, mc, stem):
    """Closed-form band (solid) against the simulated band (hatched), same
    visual language as plot_freq_term_mode_comparison's _compare figure.
    """
    z = stats.norm.ppf(1.0 - freq["AEP"].values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))

    fig, ax = plt.subplots(figsize=(11, 8.4))
    ax.fill_between(z, mc["closed_lower_cfs"], mc["closed_upper_cfs"],
                    color=C_REG, alpha=0.22, zorder=1, lw=0,
                    label="closed form (RSS of sigmas)")
    ax.fill_between(z, mc["mc_lower_cfs"], mc["mc_upper_cfs"],
                    facecolor="none", edgecolor="#1b1b1b", hatch="//",
                    lw=0.9, zorder=2,
                    label="Monte Carlo (%d draws/AEP)" % MC_N_DRAWS)
    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.2, zorder=5, label="Unregulated")
    ax.plot(z, reg_curve, color=C_REG, lw=2.6, zorder=5, label="Regulated")
    _freqmode_axes(ax)
    ax.set_title("Castle Rock regulated frequency, %d%% band\n"
                "closed-form two-piece combination vs. Monte Carlo" % pct,
                fontsize=11)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    fig.tight_layout()
    fig.savefig("%s_montecarlo_compare.png" % stem, dpi=150)
    plt.close(fig)


def write_monte_carlo_report(mc, out_path):
    """Plain-text summary: does the closed-form band hold up against a
    direct simulation of the same two sources? Written for the same
    paraphrase-into-the-memo use as freq_term_mode_comparison.txt.
    """
    lines = []

    def w(s=""):
        lines.append(s)

    worst_lo = mc["lower_diff_pct"].abs().max()
    worst_hi = mc["upper_diff_pct"].abs().max()
    worst = max(worst_lo, worst_hi)

    w("REGULATED FLOW UNCERTAINTY -- MONTE CARLO CHECK")
    w("=" * 64)
    w("Generated by #Unreg_Reg_Curve.py. Re-run the script to refresh.")
    w("")
    w("WHAT THIS CHECKS")
    w("-" * 64)
    w("The adopted band combines the frequency-term sigma and the")
    w("transform-term sigma by root-sum-of-squares on each side separately")
    w("(the closed-form two-piece-lognormal treatment). That combination")
    w("rule is exact for two ordinary Normals, but only an approximation")
    w("once the two things being combined are themselves asymmetric -- the")
    w("true combined distribution of two independent split-Normals is not")
    w("itself exactly a split-Normal. This check draws directly from the")
    w("two source distributions (%d draws per AEP, split-Normal on each"
      % MC_N_DRAWS)
    w("side using the same sigmas already adopted), pushes the frequency")
    w("draws through the actual fitted transform curve, and compares the")
    w("simulated 2.5/97.5 percentiles to the closed-form band.")
    w("")
    w("RESULT")
    w("-" * 64)
    if worst <= 5.0:
        w("The closed-form band matches the simulated band closely -- the")
        w("largest difference at any AEP is %.1f%%. The two-piece-lognormal" % worst)
        w("approximation is adequate here; no correction is warranted.")
    elif worst <= 15.0:
        w("The closed-form band is close to the simulated band but not")
        w("exact -- the largest difference at any AEP is %.1f%%, at the" % worst)
        w("most asymmetric AEPs. Small enough to note rather than act on,")
        w("but the closed-form band is an approximation, not the final word,")
        w("at the most lopsided AEPs.")
    else:
        w("The closed-form band diverges from the simulated band by more")
        w("than 15%% at at least one AEP (largest: %.1f%%). At the most" % worst)
        w("asymmetric AEPs the two-piece-lognormal approximation is being")
        w("asked to do more than it can -- the Monte Carlo result is the")
        w("more defensible number there, not the closed-form formula.")
    w("")
    w("PER-AEP COMPARISON")
    w("-" * 64)
    hdr = ("%8s  %14s  %22s  %22s  %9s  %9s" %
          ("AEP", "Reg best", "closed-form lo - hi", "Monte Carlo lo - hi",
           "% lo", "% hi"))
    w(hdr)
    w("-" * len(hdr))
    for _, r in mc.iterrows():
        w("%8.4f  %14s  %22s  %22s  %+8.1f%%  %+8.1f%%" % (
            r["AEP"], format(int(round(r["reg_best_cfs"])), ","),
            "%s - %s" % (format(int(round(r["closed_lower_cfs"])), ","),
                        format(int(round(r["closed_upper_cfs"])), ",")),
            "%s - %s" % (format(int(round(r["mc_lower_cfs"])), ","),
                        format(int(round(r["mc_upper_cfs"])), ",")),
            r["lower_diff_pct"], r["upper_diff_pct"]))
    w("")
    w("% columns = how much wider (positive) or narrower (negative) the")
    w("closed-form bound is than the Monte Carlo bound, at that AEP.")
    w("")
    w("FIGURE")
    w("-" * 64)
    w("  unreg_reg_montecarlo_compare.png  closed-form band vs. Monte Carlo")
    w("                                    band, overlaid")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def plot_final_uncertainty(freq, fit, unc, reg_curve, table_2009, stem):
    """THE adopted figure: both curves with the combined uncertainty band.

    No scatter, and no 2009 curve unless FINAL_SHOW_2009 -- everything that
    placed these lines is on unreg_reg_frequency.png. This one shows the result
    and how well it is known.
    """
    fig, ax = plt.subplots(figsize=(11, 8.4))
    z = stats.norm.ppf(1.0 - freq["AEP"].values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))

    if FINAL_SHOW_UNREG_BAND:
        ax.fill_between(z, unc["unreg_lower"], unc["unreg_upper"], color=C_UNREG,
                        alpha=0.13, zorder=1, lw=0,
                        label="Unregulated, %d%% (HEC-SSP)" % pct)
    ax.fill_between(z, unc["reg_lower"], unc["reg_upper"], color=C_REG,
                    alpha=0.17, zorder=2, lw=0,
                    label="Regulated, %d%% (frequency + transform)" % pct)
    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.6, zorder=5, label="Unregulated")
    ax.plot(z, reg_curve, color=C_REG, lw=2.6, zorder=5, label="Regulated")
    if FINAL_SHOW_2009 and table_2009 is not None and len(table_2009):
        ax.plot(stats.norm.ppf(1.0 - table_2009["AEP"].values),
                table_2009["cfs"].values, color=C_2009, lw=1.7, ls="--",
                zorder=4, label=CURVE_2009_LABEL)

    if FINAL_SHOW_SUPPORT_MARKER:
        supported = unreg_curve <= fit["x_max"]
        if supported.any() and not supported.all():
            z_edge = z[supported].max()
            ax.axvline(z_edge, color="0.35", lw=1.0, ls=":", zorder=3)
            ax.text(z_edge, FLOW_LIMITS[1] * 0.94,
                    "  transform supported to here\n  (unreg %s cfs)"
                    % format(int(fit["x_max"]), ","),
                    fontsize=8, color="0.3", va="top")

    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    # Ticks trimmed to the axis, so a label is never assigned to a tick that
    # is off the plot -- probability_axis pairs labels to ticks positionally.
    ticks = [t for t in AEP_TICKS
             if FINAL_AEP_LIMITS[1] <= t <= FINAL_AEP_LIMITS[0]]
    probability_axis(ax, ticks, FINAL_AEP_LIMITS)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title(FINAL_TITLE, fontsize=12)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    if FINAL_SHOW_FORMULA_NOTE:
        ax.text(0.995, 0.015,
                "Regulated band = RegBest +/- sqrt(freq term$^2$ + transform "
                "term$^2$), FREQ_TERM_MODE=%s\ntransform sigma hi %.3f-%.3f, "
                "lo %.3f-%.3f dex (%s, %s)"
                % (unc["freq_term_mode"],
                   np.nanmin(unc["sigma_transform_hi"]),
                   np.nanmax(unc["sigma_transform_hi"]),
                   np.nanmin(unc["sigma_transform_lo"]),
                   np.nanmax(unc["sigma_transform_lo"]),
                   TRANSFORM_SIGMA_MODE, TRANSFORM_UNCERTAINTY_BASIS),
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color="0.25",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6",
                          alpha=0.9, lw=0.8))
    fig.tight_layout()
    fig.savefig("%s_final_uncertainty.png" % stem, dpi=150)
    plt.close(fig)


def plot_final_curves(freq, fit, reg_curve, unc, table_2009, stem):
    """Regulated and unregulated lines against the 2009 study, with the
    ADOPTED uncertainty band -- the same one on unreg_reg_final_uncertainty.png.

    No scatter. Everything that went into placing these two lines is on
    unreg_reg_frequency.png; this is the one that goes in the report, where the
    points are a distraction from the result.

    reg_curve and unc are passed in rather than recomputed here, for two
    reasons found together when a user reported this figure and
    unreg_reg_final_uncertainty.png showing what looked like two different
    regulated curves:

    1. reg_curve used to be recomputed locally as apply_transform(fit,
       unreg_curve). That happened to be numerically identical to the
       reg_curve plot_frequency() returns and main() threads into
       plot_final_uncertainty() -- same fit, same unreg_curve -- but two
       independent call sites computing "the same" curve is exactly the
       kind of duplication that goes stale silently the next time only one
       of them gets edited. It is now threaded through instead.
    2. The band was the real mismatch. This figure used to draw its own
       crude +/- 1 log-sigma ribbon (transform scatter only, no frequency
       term, no z, no asymmetry) while unreg_reg_final_uncertainty.png draws
       the actual adopted combine_uncertainty() band. Two figures in one
       memo, both a red band around the same red curve, at two very
       different widths, reads as the regulated curve disagreeing with
       itself. They were never meant to be the same number -- but nothing
       distinguished them as different things, so it read as a bug. This
       figure now draws the SAME unc["reg_lower"]/unc["reg_upper"] band.
    """
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    aep_values = freq["AEP"].values
    z = stats.norm.ppf(1.0 - aep_values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))

    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.6, zorder=4,
            label="Unregulated")
    ax.plot(z, reg_curve, color=C_REG, lw=2.6, zorder=4, label="Regulated")
    ax.fill_between(z, unc["reg_lower"], unc["reg_upper"], color=C_REG,
                    alpha=0.17, zorder=1,
                    label="Regulated, %d%% (frequency + transform)" % pct)
    if table_2009 is not None and len(table_2009):
        ax.plot(stats.norm.ppf(1.0 - table_2009["AEP"].values),
                table_2009["cfs"].values, color=C_2009, lw=1.7, ls="--",
                zorder=3, label=CURVE_2009_LABEL)

    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    probability_axis(ax, AEP_TICKS, AEP_LIMITS)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("Castle Rock peak flow frequency, regulated and unregulated",
                 fontsize=12)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    fig.tight_layout()
    fig.savefig("%s_final_curves.png" % stem, dpi=150)
    plt.close(fig)


def _freqmode_bands(unc):
    """Pull out the literal and transform_curve bands from unc by name.

    combine_uncertainty() always computes both -- FREQ_TERM_MODE only picks
    which one is "reg_lower"/"reg_upper" (the adopted band) versus
    "reg_lower_alt"/"reg_upper_alt" (kept for comparison). This unwinds that
    so the comparison plots are correct no matter which mode is currently
    adopted, rather than assuming "literal" is always unc["reg_..."].
    """
    if unc["freq_term_mode"] == "literal":
        literal = (unc["reg_lower"], unc["reg_upper"])
        curve = (unc["reg_lower_alt"], unc["reg_upper_alt"])
    else:
        curve = (unc["reg_lower"], unc["reg_upper"])
        literal = (unc["reg_lower_alt"], unc["reg_upper_alt"])
    return literal, curve


def _freqmode_axes(ax):
    """Shared axis dressing for the FREQ_TERM_MODE comparison figures."""
    ax.set_yscale("log")
    ax.set_ylim(FLOW_LIMITS)
    ticks = [t for t in AEP_TICKS if FINAL_AEP_LIMITS[1] <= t <= FINAL_AEP_LIMITS[0]]
    probability_axis(ax, ticks, FINAL_AEP_LIMITS)
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: format(int(v), ",")))
    ax.grid(which="major", alpha=0.45, lw=0.8)
    ax.grid(which="minor", alpha=0.2, lw=0.5)
    ax.set_ylabel("Peak flow (cfs)")


def plot_freq_term_mode_comparison(freq, fit, reg_curve, unc, stem):
    """The two ways of combining the frequency term, plotted apart and together.

    Not a memo figure. This exists so the choice between FREQ_TERM_MODE =
    "literal" (the reviewer's formula, applied to unregulated cfs directly)
    and "transform_curve" (the same bound pushed through the fitted transform
    first, so both terms are in regulated cfs before combining) can be SEEN
    rather than argued about in the abstract. See PART 4 in the module
    docstring and write_freq_term_mode_description() for the reasoning; this
    only draws it.

    Three files:
      _freqmode_literal.png          curve + the literal band alone
      _freqmode_transform_curve.png  curve + the transform_curve band alone
      _freqmode_compare.png          both bands on one axis, one solid, one
                                      hatched, so the gap between them is
                                      visible directly rather than inferred
                                      from two separate figures
    """
    z = stats.norm.ppf(1.0 - freq["AEP"].values)
    unreg_curve = freq[FREQ_VALUE_COL].values
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))
    (literal_lo, literal_hi), (curve_lo, curve_hi) = _freqmode_bands(unc)

    panels = [
        ("literal", literal_lo, literal_hi,
         "FREQ_TERM_MODE = 'literal' -- reviewer's formula, "
         "unregulated-cfs delta added directly"),
        ("transform_curve", curve_lo, curve_hi,
         "FREQ_TERM_MODE = 'transform_curve' -- unregulated bound pushed "
         "through the fitted transform first"),
    ]
    for key, lo, hi, subtitle in panels:
        fig, ax = plt.subplots(figsize=(11, 8.4))
        ax.fill_between(z, lo, hi, color=C_REG, alpha=0.18, zorder=1, lw=0,
                        label="Regulated, %d%% band" % pct)
        ax.plot(z, unreg_curve, color=C_UNREG, lw=2.2, zorder=5,
                label="Unregulated")
        ax.plot(z, reg_curve, color=C_REG, lw=2.6, zorder=5, label="Regulated")
        _freqmode_axes(ax)
        ax.set_title("Castle Rock regulated frequency\n%s" % subtitle,
                    fontsize=11)
        ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
        fig.tight_layout()
        fig.savefig("%s_freqmode_%s.png" % (stem, key), dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 8.4))
    ax.fill_between(z, literal_lo, literal_hi, color=C_REG, alpha=0.22,
                    zorder=1, lw=0, label="literal (reviewer's formula)")
    ax.fill_between(z, curve_lo, curve_hi, facecolor="none",
                    edgecolor="#1b1b1b", hatch="//", lw=0.9, zorder=2,
                    label="transform_curve (dimensionally consistent)")
    ax.plot(z, unreg_curve, color=C_UNREG, lw=2.2, zorder=5, label="Unregulated")
    ax.plot(z, reg_curve, color=C_REG, lw=2.6, zorder=5, label="Regulated")
    _freqmode_axes(ax)
    ax.set_title("Castle Rock regulated frequency, %d%% band\n"
                "FREQ_TERM_MODE comparison: literal vs. transform_curve"
                % pct, fontsize=11)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
    fig.tight_layout()
    fig.savefig("%s_freqmode_compare.png" % stem, dpi=150)
    plt.close(fig)


def write_freq_term_mode_description(freq, unc, reg_curve, out_path):
    """Short, memo-ready methodology text -- meant to be paraphrased into
    Section 5.6, not read as an internal justification document. The longer
    reasoning (why "literal" was rejected, the full per-AEP comparison) is
    kept here too, but after the part meant for the memo, and trimmed to
    what an audit trail needs rather than the full argument that produced
    the decision.
    """
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))
    z = unc["z"]
    (literal_lo, literal_hi), (curve_lo, curve_hi) = _freqmode_bands(unc)
    aep = freq["AEP"].values

    lines = []

    def w(s=""):
        lines.append(s)

    w("REGULATED FLOW UNCERTAINTY -- METHODOLOGY")
    w("=" * 64)
    w("Generated by #Unreg_Reg_Curve.py. Re-run the script to refresh these")
    w("numbers if the analysis changes.")
    w("")
    w("FOR THE MEMO (Section 5.6)")
    w("-" * 64)
    w("The regulated flow at a given AEP carries uncertainty from two")
    w("independent sources: how well the unregulated frequency curve is")
    w("known (HEC-SSP confidence limits), and how much storm shape affects")
    w("attenuation (scatter of the regulated points about the fitted")
    w("transform). Per EM 1110-2-1619 (29 Sep 2025) Sec. 4-4.b(3), Eq. 4-6,")
    w("independent standard deviations combine by root-sum-of-squares, with")
    w("the upper and lower side of each source kept separate (Sec. 4-6a) so")
    w("the real asymmetry in both sources survives the combination.")
    w("")
    w("Before combining, the frequency term is converted into regulated-cfs")
    w("units by evaluating the fitted unregulated-to-regulated transform at")
    w("the unregulated confidence bound -- the same curve used for the best")
    w("estimate -- rather than treated as an unregulated-cfs quantity added")
    w("directly to the regulated-cfs transform term. This keeps both terms")
    w("in the same units before they are combined:")
    w("")
    w("  freq_term_hi = Transform(Unreg_97.5) - Transform(Unreg_best)")
    w("  freq_term_lo = Transform(Unreg_best) - Transform(Unreg_2.5)")
    w("  transform_term_hi/lo = the fitted curve's own confidence spread,")
    w("                         evaluated at Unreg_best")
    w("")
    w("  Upper = RegBest + sqrt(freq_term_hi^2 + transform_term_hi^2)")
    w("  Lower = RegBest - sqrt(freq_term_lo^2 + transform_term_lo^2)")
    w("")
    w("at the %d%% two-sided level (z = %.3f)." % (pct, z))
    w("")
    w("WHY NOT THE REVIEWER'S FORMULA READ LITERALLY")
    w("-" * 64)
    w("Read literally, the reviewer's formula subtracts unregulated flows")
    w("directly (Unreg_97.5 - Unreg_best) and adds that raw cfs value to the")
    w("regulated-cfs transform term -- implicitly treating a 1-cfs change in")
    w("unregulated flow as a 1-cfs change in regulated flow. That assumption")
    w("is false through most of this curve, which reduces peak flow 30-37")
    w("percent. It happens to land within a few percent of the curve-based")
    w("result over much of this particular transform's range (Appendix:")
    w("comparison table below), but that agreement is a property of this")
    w("transform's shape, not something a different basin or a re-fit")
    w("transform would be guaranteed to reproduce. Evaluating the actual")
    w("fitted curve removes the assumption instead of relying on it holding")
    w("by coincidence.")
    w("")
    w("=" * 64)
    w("AUDIT TRAIL -- not for the memo")
    w("=" * 64)
    w("")
    w("QUANTITATIVE COMPARISON")
    w("-" * 64)
    w("Adopted band (transform_curve) against the rejected literal reading,")
    w("at the ordinates the memo tabulates in Table 5-2:")
    w("")
    hdr = ("%8s  %12s  %12s  %22s  %22s  %9s" %
          ("AEP", "Unreg cfs", "Reg cfs", "transform_curve lo - hi",
           "literal lo - hi", "% diff"))
    w(hdr)
    w("-" * len(hdr))
    show_aep = [0.500, 0.100, 0.020, 0.010, 0.005, 0.002, 0.001, 0.0001]
    for a in show_aep:
        i = int(np.argmin(np.abs(aep - a)))
        u = freq[FREQ_VALUE_COL].values[i]
        r = reg_curve[i]
        ll, lh = literal_lo[i], literal_hi[i]
        cl, ch = curve_lo[i], curve_hi[i]
        pdiff = 100.0 * (lh - ch) / ch if ch else float("nan")
        w("%8.4f  %12s  %12s  %22s  %22s  %+8.1f%%" % (
            a, format(int(round(u)), ","), format(int(round(r)), ","),
            "%s - %s" % (format(int(round(cl)), ","), format(int(round(ch)), ",")),
            "%s - %s" % (format(int(round(ll)), ","), format(int(round(lh)), ",")),
            pdiff))
    w("")
    w("% diff = how much wider (positive) or narrower (negative) the")
    w("rejected literal upper bound is than the adopted transform_curve")
    w("upper bound, at that AEP.")
    w("")
    w("FIGURES")
    w("-" * 64)
    w("  unreg_reg_freqmode_literal.png          literal band alone")
    w("  unreg_reg_freqmode_transform_curve.png  transform_curve band alone")
    w("                                          (= the adopted band)")
    w("  unreg_reg_freqmode_compare.png          both bands on one axis")
    w("")
    w("FREQ_TERM_MODE is currently %r (near the top of #Unreg_Reg_Curve.py)."
      % unc["freq_term_mode"])

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


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

    synth = load_synthetic_points(SYNTH_RESULTS_CSV)
    if synth is not None:
        print("Synthetic : %d routed members, unregulated %s to %s cfs"
              % (len(synth), format(int(synth[SYNTH_UNREG_COL].min()), ","),
                 format(int(synth[SYNTH_UNREG_COL].max()), ",")))

    report_synth_by_event(synth)
    fit_x, fit_y, used = assemble_fit_points(data, wcm, synth)
    fit = build_transform(fit_x, fit_y, TRANSFORM_METHOD)

    freq = pd.read_csv(UNREG_FREQ_CSV)
    freq = freq[freq["Duration"] == FREQ_DURATION].sort_values("AEP", ascending=False)
    freq = freq.dropna(subset=["AEP", FREQ_VALUE_COL])
    table_2009 = curve_2009_frame()

    print("=" * 78)
    print("Transform : %s   n = %d" % (TRANSFORM_METHOD, fit["n"]))
    print("Line drawn from:")
    for source in ("adjusted", "wcm_rc", "synthetic"):
        on = FIT_SOURCES.get(source, False)
        note = ("%d points" % used[source]) if source in used else (
            "on, but no points available" if on else "shown only, does not move the line")
        print("             [%s] %-10s %s" % ("x" if on else " ", source, note))
    if not FIT_SOURCES.get("wcm_rc", False):
        print("             wcm_rc is the pre-adjustment set -- no Obs_RC data")
        print("             exists behind it, so it carries an uncorrected bias.")
    if TRANSFORM_METHOD == "loess":
        print("            LOESS span %.2f, pseudo r2 = %.4f" % (LOESS_SPAN, fit["r2"]))
        print("            power law for reference: reg = %.4g x unreg^%.4f, r2 = %.4f"
              % (fit["power"]["a"], fit["power"]["b"], fit["power"]["r2"]))
    else:
        print("            reg = %.4g x unreg^%.4f   log r2 = %.4f"
              % (fit["power"]["a"], fit["power"]["b"], fit["power"]["r2"]))
    print("Supported over unregulated %s to %s cfs"
          % (format(int(fit["x_min"]), ","), format(int(fit["x_max"]), ",")))
    if TRANSFORM_SIGMA_MODE == "local":
        edges = transform_sigma_dex(fit, np.array([fit["x_min"], fit["x_max"]]))
        print("Scatter about the line: x/ %.3f at %s cfs to x/ %.3f at %s cfs"
              % (10 ** edges[0], format(int(fit["x_min"]), ","),
                 10 ** edges[1], format(int(fit["x_max"]), ",")))
        print("            (local, span %.2f; pooled would be x/ %.3f everywhere)"
              % (TRANSFORM_SIGMA_SPAN, 10 ** fit["se_dex"]))
    else:
        print("Scatter about the line: x/ %.3f (1 sigma, pooled)"
              % 10 ** fit["se_dex"])
    print("Monotonic enforced: %s   clipped at 1:1: %s"
          % (ENFORCE_MONOTONIC, CLIP_TO_UNREG))
    print("=" * 78)

    plot_scatter(data, fit, wcm, synth, PLOT_STEM)
    reg_curve = plot_frequency(freq, data, fit, wcm, synth, table_2009,
                               PLOT_STEM)

    # unc is computed before plot_final_curves (moved up from after it) so
    # that figure can draw the SAME adopted band as plot_final_uncertainty,
    # instead of the crude, unrelated "+/- 1 std error" ribbon it used to
    # draw locally. Two figures in the same memo both showing a red band
    # around the regulated curve, at two different widths, read as the
    # regulated curve itself disagreeing between figures -- which is what
    # got reported. They were never the same quantity: one was the full
    # combined uncertainty band, the other a rough single-sigma diagnostic
    # nobody meant to put next to the adopted band. Now there is one band.
    unc = combine_uncertainty(freq, fit, reg_curve)
    report_uncertainty(freq, unc, fit)

    if MAKE_FINAL_CURVE_PLOT:
        plot_final_curves(freq, fit, reg_curve, unc, table_2009, PLOT_STEM)
    if MAKE_FINAL_UNCERTAINTY_PLOT:
        plot_final_uncertainty(freq, fit, unc, reg_curve, table_2009, PLOT_STEM)
    if MAKE_FREQ_MODE_COMPARISON:
        plot_freq_term_mode_comparison(freq, fit, reg_curve, unc, PLOT_STEM)
        desc_path = os.path.join(os.path.dirname(PLOT_STEM),
                                 "freq_term_mode_comparison.txt")
        write_freq_term_mode_description(freq, unc, reg_curve, desc_path)
        print("FreqMode: %s_freqmode_literal.png, %s_freqmode_transform_curve.png,"
              % (PLOT_STEM, PLOT_STEM))
        print("          %s_freqmode_compare.png" % PLOT_STEM)
        print("          %s" % desc_path)

    if RUN_MONTE_CARLO_CHECK:
        mc = monte_carlo_check(freq, fit, unc, reg_curve)
        plot_monte_carlo_check(freq, fit, reg_curve, mc, PLOT_STEM)
        mc_path = os.path.join(os.path.dirname(PLOT_STEM), "monte_carlo_check.txt")
        write_monte_carlo_report(mc, mc_path)
        worst = max(mc["lower_diff_pct"].abs().max(), mc["upper_diff_pct"].abs().max())
        print("MonteCarlo: %s_montecarlo_compare.png" % PLOT_STEM)
        print("            %s" % mc_path)
        print("            largest closed-form vs. MC difference: %.1f%%" % worst)

    out = freq[["AEP", "Value", FREQ_VALUE_COL]].copy()
    out = out.rename(columns={"Value": "unreg_computed_cfs",
                              FREQ_VALUE_COL: "unreg_expected_cfs"})
    band = 10 ** transform_sigma_dex(fit, out["unreg_expected_cfs"].values)
    out["reg_inferred_cfs"] = reg_curve
    out["reg_lower_1se_cfs"] = reg_curve / band
    out["reg_upper_1se_cfs"] = reg_curve * band
    # --- combined uncertainty (PART 4), reviewer's formula ------------------
    pct = int(round(100 * UNCERTAINTY_CONF_LEVEL))
    out["unreg_lower_%dpct_cfs" % pct] = unc["unreg_lower"]
    out["unreg_upper_%dpct_cfs" % pct] = unc["unreg_upper"]
    out["reg_lower_%dpct_cfs" % pct] = unc["reg_lower"]
    out["reg_upper_%dpct_cfs" % pct] = unc["reg_upper"]
    out["freq_term_mode"] = unc["freq_term_mode"]
    out["transform_slope_b"] = unc["slope"]
    out["sigma_freq_lo_dex"] = unc["sigma_freq_lo"]
    out["sigma_freq_hi_dex"] = unc["sigma_freq_hi"]
    out["sigma_transform_lo_dex"] = unc["sigma_transform_lo"]
    out["sigma_transform_hi_dex"] = unc["sigma_transform_hi"]
    out["freq_term_lo_cfs"] = unc["freq_term_lo_cfs"]
    out["freq_term_hi_cfs"] = unc["freq_term_hi_cfs"]
    out["transform_term_lo_cfs"] = unc["transform_term_lo_cfs"]
    out["transform_term_hi_cfs"] = unc["transform_term_hi_cfs"]
    out["delta_lo_cfs"] = unc["delta_lo_cfs"]
    out["delta_hi_cfs"] = unc["delta_hi_cfs"]
    out["band_asymmetry"] = unc["asymmetry"]
    out["reg_upper_clipped_at_unreg"] = unc["reg_upper_clipped"]
    # --- the OTHER FREQ_TERM_MODE, for comparison only -----------------------
    out["alt_freq_term_mode"] = unc["alt_freq_term_mode"]
    out["reg_lower_%dpct_alt_cfs" % pct] = unc["reg_lower_alt"]
    out["reg_upper_%dpct_alt_cfs" % pct] = unc["reg_upper_alt"]
    out["freq_term_lo_literal_cfs"] = unc["freq_term_lo_literal_cfs"]
    out["freq_term_hi_literal_cfs"] = unc["freq_term_hi_literal_cfs"]
    out["freq_term_lo_curve_cfs"] = unc["freq_term_lo_curve_cfs"]
    out["freq_term_hi_curve_cfs"] = unc["freq_term_hi_curve_cfs"]
    out["reg_powerlaw_cfs"] = apply_power_law(fit["power"],
                                              out["unreg_expected_cfs"].values)
    out["reduction_pct"] = 100.0 * (1.0 - reg_curve / out["unreg_expected_cfs"])
    out["extrapolated"] = out["unreg_expected_cfs"] > fit["x_max"]
    out["reg_2009_cfs"] = interp_2009(out["AEP"].values, table_2009)
    out["reg_minus_2009_cfs"] = out["reg_inferred_cfs"] - out["reg_2009_cfs"]
    out["reg_vs_2009_pct"] = 100.0 * (out["reg_inferred_cfs"] / out["reg_2009_cfs"]
                                      - 1.0)
    out.to_csv(os.path.join(OUT_DIR, "regulated_frequency_inferred.csv"),
               index=False, float_format="%.8g")
    # NOT "%.1f". One format has to serve flows in the hundreds of thousands
    # and log-space sigmas of a few hundredths, and a single decimal place
    # destroys the small ones: it wrote every AEP below 0.05 as 0.0, turned
    # 0.99 into 1.0 and 0.95 into 0.9, and flattened every sigma and the
    # transform slope to 0.0 or 0.1. Significant figures, not decimal places.

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
    if MAKE_FINAL_CURVE_PLOT:
        print("Final : %s_final_curves.png  (lines only, no points)" % PLOT_STEM)
    if MAKE_FINAL_UNCERTAINTY_PLOT:
        print("ADOPTED: %s_final_uncertainty.png  (both curves + %d%% band%s)"
              % (PLOT_STEM, int(round(100 * UNCERTAINTY_CONF_LEVEL)),
                 "" if FINAL_SHOW_2009 else ", 2009 curve off"))
    print("Table : %s/regulated_frequency_inferred.csv" % OUT_DIR)


main()
