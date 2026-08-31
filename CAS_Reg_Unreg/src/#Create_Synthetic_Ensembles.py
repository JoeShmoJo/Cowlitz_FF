#Create_Synthetic_Ensembles.py
# -*- coding: utf-8 -*-
"""
Build scaled synthetic flood ensembles for the unregulated-regulated transform.

WHY
---
Only one water year in 99 exceeds the unregulated 100-year peak, so the top of
the unreg-reg curve is essentially unconstrained by the record. These synthetics
populate that range by scaling observed events up to target magnitudes.

DESIGN -- a full factorial over the three things that could matter:

    SHAPE     4 source events spanning sharp to sustained hydrographs.
              Attenuation among large observed events is NOT predictable from
              shape, magnitude, antecedent flow or timing (all r-squared < 0.10,
              n=18), so shape is varied rather than assumed away.

    MAGNITUDE 4 targets: the 100-year, 250-year and 500-year unregulated peaks,
              plus one BEYOND the 500-year. The last one exists so the drawn
              curve is not extrapolating past its final point exactly where the
              answer is needed.

    POOL      starting pool elevation. FIVE bases are available and each is
              switched on or off individually in POOL_BASES_ENABLED: the WCM
              rule curve, the 50% pool duration curve for that calendar date,
              the observed pool on that date, the median pool from the
              unregulated POR run, and the highest of whichever are on. Each
              enabled basis is written as its own member so the results can be
              compared, or pooled into a Monte Carlo over starting conditions.

    -> the member count is 4 events x 4 magnitudes x (bases switched on).
       With only the rule curve on, that is 16 members; with all five, 80.
       The script prints the on/off list and the resulting count at the top of
       its run, so what was actually built is never in doubt.

       A basis that is enabled but has no value for an event (no POR run, no
       observed record that far back) falls back to POOL_FALLBACK_BASIS, and
       the member is tagged in pool_basis_used rather than substituted quietly.

SCALING
-------
BOTH the Mossyrock inflow and the Castle Rock local are multiplied by the SAME
factor. The local is a median 48% of unregulated flow at large events, and the
project has no control over it, so scaling only the reservoir inflow would
shrink the uncontrolled fraction and change attenuation for reasons unrelated to
magnitude. Uniform scaling preserves the observed coincidence between the two.

The scale factor is set so the scaled unregulated PEAK (Mossyrock inflow plus
local, the same quantity the transform uses on its x-axis) hits the target.

Outputs: the ensemble DSS, a mapping CSV carrying scale factor / target AEP /
source event / pool basis for every member, and a diagnostic plot per event.

!! THE OUTPUT DSS IS HAND-EDITED AFTER THIS SCRIPT RUNS !!
---------------------------------------------------------
Matching a peak AND a 5-day volume on a sharp source event forces f_out above
f_peak (1.9x vs 1.1x for Dec1977), so a small rise the observed hydrograph
already had gets stretched into what looks like a second flood a day or two
ahead of the main peak. On the Dec1977 and Nov1986 members it reached 0.43-0.45
of the peak, and both are CHOPPED BY HAND IN DSSVUE -- on the MOSSYROCK/FLOW-IN
and CASTLE ROCK/FLOW-LOCAL records -- before the ensemble is loaded into ResSim.
See "The hand-chopped synthetics" in RUN_ORDER.md.

build_set DELETES OUT_DSS before writing, so RUNNING THIS SCRIPT DESTROYS THAT
EDIT. Nothing downstream catches it: the members still hit their peak and volume
targets with the bump present, so the next ResSim run just quietly routes the
two-peaked hydrographs again. If this script is re-run, the chop has to be
redone before step 23.

Damping the bump in code was tried in Aug 2026 and reverted -- the safe amount
of damping (past ~1.4 a shoulder scales above the peak and invents a new
maximum) only moved it 1-2%.
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
OBS_DSS = r"../../CAS_Unreg_FF/data/obsData.dss"

# EXTERNAL: requires the ResSim watershed (not in this repository)
OUT_DSS = (r"../output/ensemble_synthetic.dss")
OUT_DSS_VERSION = 7

OUT_DSS_ALT = OUT_DSS.replace(".dss", "_lineartaper.dss")
MAPPING_CSV = r"../output/ensemble_synthetic_mapping.csv"
PLOT_STEM = r"../output/diagnostics/ensemble_synthetic"

PATH_MOS_IN = "//MOSSYROCK/FLOW-IN/*/1HOUR/FOR_RESSIM/"
PATH_CAS_LOCAL = "//CASTLE ROCK/FLOW-LOCAL/*/1HOUR/FOR_RESSIM/"
PATH_MOS_ELEV_DAILY = "//MOS/ELEV/*/1DAY/USGS/"

ENS_SUFFIX = "SYNTH"

# --- second plot set: the inflows ResSim ACTUALLY ran ------------------------
# The _events.png plot recomputes each member from the scaling functions, so it
# shows the ensemble AS BUILT. That is not what was routed: the Dec1977 and
# Nov1986 mini peaks are chopped by hand in DSSVue after this script writes the
# DSS (see the warning at the top of this file), and the plot has no way to
# know. The figure therefore shows two-peaked hydrographs that no longer exist
# in the run the results came from.
#
# When SIM_DSS is present, a second set is drawn straight from the simulation's
# own inflow records, which are the hand-edited ones ResSim read. That set is
# the one to put in a memo.
#
# The observed source hydrograph is not in the simulation -- ResSim never saw
# the unscaled event -- so it still comes from IN_DSS, the same place the build
# takes it from.
SIM_DSS = r"../output/ResSim_Synth.dss"
PLOT_AS_RUN = True
# Pathname parts in the simulation output. ResSim writes its own B parts, which
# are not the ones the ensemble was written under.
SIM_PART_B_MOS = "MOSSYROCK-POOL"
SIM_PART_C_MOS = "FLOW-IN"
SIM_PART_B_CAS = "CASTLEROCK_NWS"
SIM_PART_C_CAS = "FLOW-LOCAL"
SIM_PART_F = "RESSIM_SYNTH"
# Overlay the as-built curve as a thin dashed line, so the hand edit shows up
# as the difference between the two. Useful for checking the chop; noise on a
# figure meant for a memo, which is why it is off.
ASRUN_SHOW_BUILT = False
# A member whose largest hourly difference from the as-built curve exceeds this
# percentage of its own peak is reported as edited. Small non-zero differences
# are normal -- the simulation writes what it read at its own precision.
EDIT_DETECT_PCT = 2.0
# Hours to shift the simulation record onto the build index. None = detect it
# from the data (see detect_sim_hour_offset); an integer forces it.
SIM_HOUR_OFFSET = None
SIM_OFFSET_SEARCH = 4

# Rebuild the plots from the existing mapping CSV and the simulation WITHOUT
# touching the ensemble DSS.
#
# Turn this on to regenerate a figure. build_set DELETES OUT_DSS before writing,
# so a normal run to refresh a plot would destroy the hand-chopped members --
# the exact trap described at the top of this file. With PLOTS_ONLY the DSS is
# never opened, and the mapping is read back from MAPPING_CSV instead of being
# rebuilt.
PLOTS_ONLY = True

# DQC: the 12-panel figure is not legible at memo size. Two summary figures
# are drawn alongside it -- every source storm overlaid unscaled, and the full
# scaling family for one representative member. The 12-panel set moves to an
# appendix at full page size.
PLOT_SUMMARY_FIGURES = True
FEATURED_EVENT = "Feb1996"      # falls back to the first enabled event

# --- source events: (label, peak date, shape 5day/peak, note) ---------------
# Window is centred on the peak; WINDOW_BEFORE/AFTER days on each side.
# --- SOURCE EVENTS -- TOGGLE HERE --------------------------------------------
# The catalog below is every water year in the record ranked by UNREGULATED peak
# at Castle Rock, from the WCM_RC period-of-record run. To change which storms
# the synthetics are built from, flip the first column. Nothing else needs
# editing, and events can be added freely -- the simulation is cheap, and more
# source events is the only way to learn whether the spread at the top of the
# regulated curve is real or an artefact of a thin sample.
#
# WHY BIG EVENTS ONLY
# Hydrograph shape is not independent of size. A large storm is large because it
# was widespread and sustained, and that is exactly what determines whether the
# reservoir can hold it. Scaling a small event up to a 500-year magnitude keeps
# the small event's shape and asks the reservoir to route a flood volume that
# storm never produced -- the routing is arithmetically fine and physically
# meaningless. The scale factor column below is the honest measure of how far
# each event is being stretched; anything much past 2x is a stretch in both
# senses. SCALE_FACTOR_WARN flags them at run time.
#
# The ratio column is the WCM_RC regulated/unregulated peak ratio -- how much
# the project actually took off that storm. It is the best available proxy for
# shape, and a good source set spans it rather than clustering: 1.00 means the
# flood passed straight through, 0.5 means half of it was held back.
#
#     use    label      peak date     unreg peak  ratio   note
SOURCE_EVENT_CATALOG = [
    (True,  "Feb1996", "1996-02-09",     212246, 0.62, "largest on record"),
    (True,  "Nov2006", "2006-11-07",     155018, 0.61, "rank 2"),
    (True,  "Dec2015", "2015-12-09",     152270, 0.53, "rank 3, well attenuated"),
    (True,  "Dec1933", "1933-12-23",     149770, 1.00, "PASS-THROUGH -- storage exhausted"),
    (True,  "Dec1977", "1977-12-02",     148667, 0.54, "rank 5"),
    (True,  "Jan2009", "2009-01-08",     143468, 0.57, "rank 6"),
    (True,  "Dec1975", "1975-12-04",     136059, 0.56, "rank 7"),
    (True,  "Jan1990", "1990-01-09",     130043, 0.65, "rank 8, mid attenuation"),
    (True, "Feb2003", "2003-02-01",     119348, 0.52, "rank 9"),
    (True, "Dec2025", "2025-12-11",     116111, 0.52, "rank 10, most recent large event"),
    (True, "Nov1990", "1990-11-25",     108150, 0.58, "rank 11"),
    (True, "Nov1986", "1986-11-24",     104890, 0.70, "rank 12, sustained"),
    (False, "Feb1986", "1986-02-24",      94558, 0.66, "rank 13"),
    (False, "Dec1996", "1996-12-30",      90906, 0.60, "rank 14"),
    (False, "Jan1972", "1972-01-20",      90280, 0.71, "rank 15, sustained"),
    (False, "Dec2007", "2007-12-04",      89749, 0.68, "rank 16"),
    (False, "Dec1946", "1946-12-13",      88542, 0.61, "rank 17"),
    (False, "Dec1964", "1964-12-23",      87903, 0.59, "rank 18"),
    (False, "Jan1982", "1982-01-24",      87607, 0.74, "rank 19, sustained"),
    (False, "Dec2023", "2023-12-06",      87399, 0.60, "rank 20"),
]
# Warn when an enabled event has to be stretched further than this to reach the
# largest target. Not an error -- a flag that the shape being routed may not
# belong to a flood of that size.
SCALE_FACTOR_WARN = 2.0
# Refuse to enable an event whose unregulated peak is below this. Set to 0 to
# allow anything.
MIN_SOURCE_PEAK_CFS = 100000.0

SOURCE_EVENTS = [(label, date, "%s cfs unreg, WCM_RC ratio %.2f -- %s"
                  % (format(peak, ","), ratio, note))
                 for use, label, date, peak, ratio, note in SOURCE_EVENT_CATALOG
                 if use]
WINDOW_BEFORE_DAYS = 10       # lead-in, so the pool responds before the peak
WINDOW_AFTER_DAYS = 20        # recession, long enough to pass the flood

# --- how the events are scaled ----------------------------------------------
# "volume_matched": the multiplier is a function of FLOW, not of time --
#     f(Q) = f_out + (f_peak - f_out) * (Q - Qmin) / (Qmax - Qmin)
#     so the peak scales by f_peak exactly and every lower flow scales by less
#     (or more) in a way that is monotone in Q. That matters: a time-based taper
#     can scale a shoulder harder than the peak and lift it ABOVE the peak,
#     inventing a new maximum and breaking the target. A flow-based multiplier
#     cannot reorder the hydrograph, so the peak stays the peak and the result
#     still looks like a flood.
#     f_out is solved so the 5-day volume about the peak hits its own target,
#     then refined by iteration because changing the shape moves which 5-day
#     window is the maximum one.
# "linear_taper": the simple fallback. The multiplier is f_peak at the peak and
#     falls linearly to 1.0 at +/- TAPER_HALF_WIDTH_DAYS, leaving everything
#     outside untouched. Only the peak target is met; volume lands where it
#     lands. Written to its own DSS so the two can be compared.
SCALING_METHOD = "volume_matched"
ALSO_WRITE_LINEAR_TAPER = True          # second DSS using the fallback method
VOLUME_HALF_WIDTH_DAYS = 2.5            # +/- window the 5-day volume is matched over
TAPER_HALF_WIDTH_DAYS = 2.5             # fallback taper half width
VOLUME_ITERATIONS = 12                  # refine f_out against the true rolling max
VOLUME_TOLERANCE = 0.002                # stop when the 5-day is within 0.2%
SHAPE_STRAIN_WARN = 2.0                 # warn when f_out/f_peak is beyond this or its inverse

# --- how fast the hydrograph returns to the observed one outside the window --
# The flow-based multiplier is a function of Q alone, so on its own it rescales
# the ENTIRE member window -- the 10-day lead-in and the 20-day recession
# included. Nothing outside +/- VOLUME_HALF_WIDTH_DAYS is constrained by either
# target, so that part of the hydrograph is being changed for no reason, and it
# changes the antecedent condition the reservoir starts the flood with.
#
# The multiplier is therefore wrapped in a time envelope:
#
#     m(Q, t) = 1 + (f(Q) - 1) * e(t)
#
#     e(t) = 1                     inside  +/- VOLUME_HALF_WIDTH_DAYS
#     e(t) -> 0 over RETURN_DAYS   outside, measured from the WINDOW EDGE
#     e(t) = 0                     beyond that -- the observed flow, untouched
#
# Inside the window nothing changes, so the peak and 5-day targets are still met
# exactly. The envelope is inside the iteration loop, so the 5-day solve sees
# the hydrograph that is actually written.
#
# RETURN_DAYS is the knob to turn. Smaller = back on the observed hydrograph
# sooner, at the cost of a sharper join at the window edge; larger = a gentler
# join that touches more of the recession. 1.5 days puts the member back on the
# observed flow four days out from the peak.
OUTSIDE_RETURN_DAYS = 1.5
# Recession side only. None = same as OUTSIDE_RETURN_DAYS. Worth lengthening if
# a step is visible on the falling limb, which is the side that carries volume.
OUTSIDE_RETURN_DAYS_AFTER = None
# Shape of the return:
#   "cosine" : half-cosine, flat at both ends -- no kink where it joins the
#              window and none where it reaches the observed flow. The default.
#   "linear" : straight ramp. Fastest to explain, slight kink at each end.
#   "power"  : (1 - d/L) ** OUTSIDE_RETURN_POWER. Exponent above 1 holds the
#              scaling longer then drops; below 1 drops immediately then trails.
OUTSIDE_RETURN_SHAPE = "cosine"
OUTSIDE_RETURN_POWER = 2.0
# False restores the old behaviour: the multiplier applies across the whole
# member window with no return at all.
APPLY_OUTSIDE_RETURN = True

# --- starting pool bases -- TOGGLE HERE --------------------------------------
# One member is written per source event x magnitude x ENABLED basis, so this
# block is what sets the member count. Turning one off removes it entirely; it
# is not written and not plotted.
#
#   "rulecurve"  the WCM rule curve on the event's start date. Always available.
#   "duration50" the 50% exceedance pool for that calendar date, from the
#                OBSERVED daily elevation record (PATH_MOS_ELEV_DAILY in
#                obsData.dss). Needs DURATION_50_MIN_YEARS of record behind that
#                calendar day.
#   "observed"   the observed pool on that calendar date, from the same record.
#   "median_por" the same statistic -- the median pool for that calendar day --
#                but computed from the unregulated period-of-record ResSim run
#                instead. Needs POR_ELEV_DSS, which is external to this
#                repository.
#   "highest"    the highest of whichever of the above are enabled AND
#                available -- the conservative convention from the workflow
#                email. Not a record in its own right; it picks one per event.
#
# NOTE: 50% exceedance, 50% non-exceedance and the median are the same number,
# so "duration50" and "median_por" are the SAME STATISTIC from two different
# sources -- observed record vs POR simulation. Switch both on to compare them;
# do not treat them as two different concepts.
#
# A basis that is enabled but unavailable for an event falls back to the rule
# curve, and the member is tagged in pool_basis_used so the substitution is on
# the record rather than silent.
POOL_BASES_ENABLED = {
    "rulecurve": False,
    "duration50": False,
    "observed": False,
    "median_por": True,
    "highest": False,
}
# Where an enabled basis has no value for an event.
POOL_FALLBACK_BASIS = "rulecurve"

# What the ELEV lookback record written into the ensemble looks like, for the
# statistic-based bases (duration50, median_por).
#   "trace" : the day-varying statistic itself -- a member whose window starts
#             15 Dec gets 15 Dec's median at its first stamp, 16 Dec's the next
#             day, and so on. The record IS the median-by-calendar-day curve.
#   "flat"  : the start date's value held constant across the whole window.
# ResSim reads only the value at the simulation start, and that value is
# IDENTICAL either way, so this does not change a run. It changes what the
# record looks like in DSSVue and what anything reading the full series sees.
# "trace" is the honest one: the record then says what it claims to be.
POOL_SERIES_STYLE = "trace"

# The order members are written in, for any basis switched on above.
POOL_BASIS_ORDER = ["rulecurve", "duration50", "observed", "median_por", "highest"]

POR_ELEV_DSS = (r"C:\Projects\2026_Cowlitz_Flow_Frequency\ResSim\NWP_CowlitzLewis"
                r"\watershed\NWP_CowlitzLewis_ResSim4\rss\Unreg_POR_FIS\simulation.dss")
POR_ELEV_PATH = "//Mossyrock-Pool/Elev/*/1Hour/Ensemble--0/"
POR_ELEV_MIN_YEARS = 10

# --- target unregulated PEAK magnitudes -------------------------------------
# From CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv, Duration = Peak.
# Peak AND 5-day targets are read from the unregulated frequency table, so the
# two constraints always come from the same curve and column.
UNREG_FREQ_CSV = r"../../CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv"
FREQ_VALUE_COL = "Expected"
TARGET_AEPS = [("100yr", 0.010), ("250yr", 0.004), ("500yr", 0.002)]
# One target beyond the reported curve, to anchor the end slope of the drawn
# curve. Given as a multiple of the smallest-AEP value on the table.
BEYOND_LABEL = "beyond"
BEYOND_FACTOR = 1.20

RULE_CURVE_ANCHORS = [
    (1, 1, 745.5), (1, 31, 745.5), (6, 1, 778.5), (9, 30, 778.5), (12, 1, 745.5),
]
LOOKBACK_DAYS = 1
ELEV_DAILY_TO_HOURLY = "interpolate"
DURATION_50_MIN_YEARS = 10     # min years of record behind a duration-curve value

# Synthetic calendar the members are labelled on (any non-leap year)
ENS_LABEL_START = datetime(1999, 10, 1, 0, 0)

# Every member shares its source event's real dates, so the 12 members built
# from one event would collide when results are reassembled. Each member is
# therefore given its own synthetic water year: the month and day of the source
# event are kept, so seasonality stays visible, and only the year is replaced.
# Years start well before the record so they can never be mistaken for real
# ones. The mapping CSV carries the synthetic date AND the true source date.
SYNTH_YEAR_BASE = 1800

CLIP_NEGATIVE_FLOW = False
SENTINEL = -901.0

# ----------------------------------------------------------------------------


def check_synthetic_calendar(mapping):
    """Verify the promise the synthetic calendar makes, instead of asserting it.

    One member per synthetic water year, each member's whole window inside that
    one water year, and no two windows sharing an hour. A failure here is what
    produces water years with two members' worth of hours, water years with
    none, and 'overlapping timestamps' warnings in the extract step.
    """
    rows = mapping.sort_values("real_start").reset_index(drop=True)
    problems = []

    straddle = rows[rows["real_start"].map(water_year)
                    != rows["real_end"].map(water_year)]
    for _, r in straddle.iterrows():
        problems.append("member %d crosses 01 Oct: %s -> %s"
                        % (r["member"], r["real_start"].date(),
                           r["real_end"].date()))

    actual_wy = rows["real_start"].map(water_year)
    mismatch = rows[actual_wy != rows["synth_water_year"]]
    for _, r in mismatch.iterrows():
        problems.append("member %d labelled WY%d but starts in WY%d (%s)"
                        % (r["member"], r["synth_water_year"],
                           water_year(r["real_start"]), r["real_start"].date()))

    for wy, grp in rows.groupby(actual_wy):
        if len(grp) > 1:
            problems.append("WY%d holds %d members: %s"
                            % (wy, len(grp), list(grp["member"])))

    for i in range(1, len(rows)):
        if rows.loc[i, "real_start"] <= rows.loc[i - 1, "real_end"]:
            overlap = (rows.loc[i - 1, "real_end"] - rows.loc[i, "real_start"])
            problems.append("members %d and %d overlap by %d hours"
                            % (rows.loc[i - 1, "member"], rows.loc[i, "member"],
                               overlap.total_seconds() / 3600 + 1))

    if problems:
        print("   *** SYNTHETIC CALENDAR PROBLEMS (%d) ***" % len(problems))
        for line in problems[:20]:
            print("      %s" % line)
        if len(problems) > 20:
            print("      ...and %d more" % (len(problems) - 20))
    else:
        print("   calendar check: %d members, %d water years, no straddles, "
              "no overlaps" % (len(rows), rows["synth_water_year"].nunique()))
    return problems


def water_year(stamp):
    """Water year containing a timestamp: 01 Oct YYYY-1 through 30 Sep YYYY."""
    stamp = pd.Timestamp(stamp)
    return stamp.year + (1 if stamp.month >= 10 else 0)


def dss_version(path):
    """DSS file version from the header: byte 12 is 6 for v6, 0 for v7."""
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
    anchored = daily.copy()
    anchored.index = anchored.index + pd.Timedelta(hours=12)
    return (anchored.reindex(anchored.index.union(hourly_index))
            .interpolate(method="time").reindex(hourly_index))


def rule_curve_on_index(index):
    """Seasonal WCM rule curve, linear between anchors.

    Interpolate on float HOURS from a fixed epoch. Do NOT mix Index.asi8 with
    Timestamp.value -- pandas 3 returns microseconds from asi8 and nanoseconds
    from .value, and the 1000x mismatch silently flattens the curve.
    """
    epoch = pd.Timestamp(1900, 1, 1)
    knots = []
    for year in range(int(index.year.min()) - 1, int(index.year.max()) + 2):
        for month, day, elev in RULE_CURVE_ANCHORS:
            knots.append(((pd.Timestamp(year, month, day) - epoch)
                          / pd.Timedelta(hours=1), float(elev)))
    knots.sort()
    kt = np.array([k[0] for k in knots], dtype=float)
    kv = np.array([k[1] for k in knots], dtype=float)
    x = (index - epoch) / pd.Timedelta(hours=1)
    return np.interp(np.asarray(x, dtype=float), kt, kv)


def duration_50_by_dayofyear(daily_elev, min_years):
    """50% exceedance pool elevation for each (month, day) from the record."""
    frame = pd.DataFrame({"elev": daily_elev.values,
                          "month": daily_elev.index.month,
                          "day": daily_elev.index.day})
    grouped = frame.groupby(["month", "day"])["elev"]
    median = grouped.median()
    counts = grouped.count()
    return median[counts >= min_years]


def dayofyear_on_index(table, index):
    """Day-of-year statistic as a series on a DatetimeIndex (Feb 29 -> Feb 28).

    Serves both the observed duration curve and the POR median, which are the
    same kind of table: keyed on (month, day).
    """
    out = np.full(len(index), np.nan)
    for i in range(len(index)):
        key = (int(index.month[i]), int(index.day[i]))
        if key in table.index:
            out[i] = table.loc[key]
        elif key == (2, 29) and (2, 28) in table.index:
            out[i] = table.loc[(2, 28)]
    return pd.Series(out, index=index).interpolate(limit_direction="both").values


def pool_series_for(table, index, start, style):
    """The ELEV lookback record for a statistic-based basis.

    "trace" writes the day-varying statistic, so the record starting on the
    event's start date carries that date's value at its first stamp and follows
    the calendar-day curve from there. "flat" holds the start value across the
    window. The value AT THE START -- the only one ResSim reads -- is the same
    either way.
    """
    start_value = pool_on_date(table, start)
    if style == "trace":
        return pd.Series(dayofyear_on_index(table, index), index=index), start_value
    return pd.Series(start_value, index=index), start_value


def read_targets(csv_path, value_col):
    """Peak and 5-day unregulated targets per AEP, from one frequency table."""
    table = pd.read_csv(csv_path)
    out = {}
    for label, aep in TARGET_AEPS:
        row = {}
        for duration, key in [("Peak", "peak"), ("5-Day", "vol5")]:
            sub = table[table["Duration"] == duration].dropna(subset=[value_col])
            pick = sub.iloc[(sub["AEP"] - aep).abs().argsort()[:1]]
            row[key] = float(pick[value_col].iloc[0])
            row["aep"] = float(pick["AEP"].iloc[0])
        out[label] = row
    smallest = min(a for _, a in TARGET_AEPS)
    base = out[[l for l, a in TARGET_AEPS if a == smallest][0]]
    out[BEYOND_LABEL] = {"peak": base["peak"] * BEYOND_FACTOR,
                         "vol5": base["vol5"] * BEYOND_FACTOR, "aep": np.nan}
    return out


def rolling_max_mean(values, index, hours):
    """Largest N-hour mean of a series (the same statistic the targets are)."""
    return float(pd.Series(values, index=index)
                 .rolling(hours, min_periods=int(hours * 0.9)).mean().max())


def return_envelope(index, peak_time, half_width_days):
    """Weight on the scaling: 1 inside the window, falling to 0 outside it.

    Distance is measured from the EDGE of the window, not from the peak, so
    half_width_days is untouched no matter what the return settings are.
    """
    offset = (index - peak_time) / pd.Timedelta(days=1)
    if not APPLY_OUTSIDE_RETURN:
        return np.ones(len(offset))
    beyond = np.abs(offset) - half_width_days
    length = np.where(offset < 0, OUTSIDE_RETURN_DAYS,
                      OUTSIDE_RETURN_DAYS if OUTSIDE_RETURN_DAYS_AFTER is None
                      else OUTSIDE_RETURN_DAYS_AFTER)
    length = np.maximum(np.asarray(length, dtype=float), 1e-9)
    d = np.clip(beyond / length, 0.0, 1.0)          # 0 at the edge, 1 when done
    if OUTSIDE_RETURN_SHAPE == "linear":
        return 1.0 - d
    if OUTSIDE_RETURN_SHAPE == "power":
        return (1.0 - d) ** OUTSIDE_RETURN_POWER
    return 0.5 * (1.0 + np.cos(np.pi * d))          # cosine, flat at both ends


def apply_multiplier(flow, factor, envelope):
    """Blend a flow-based multiplier back to 1.0 through the time envelope."""
    return (1.0 + (factor - 1.0) * envelope) * flow


def peak_correction(factor_peak, ratio):
    """Rescale that fixes the peak without disturbing the returned-to region.

    A flat multiply would lift the whole member, including the part that is
    supposed to be the observed hydrograph exactly. Correcting the MULTIPLIER
    instead leaves m = 1 wherever the envelope has already reached zero.
    """
    if not np.isfinite(ratio) or abs(factor_peak - 1.0) < 1e-9:
        return 1.0
    return (factor_peak * ratio - 1.0) / (factor_peak - 1.0)


def outside_change(scaled, flow, index, peak_time, half_width_days):
    """How much of the member outside the volume window is not the observed flow."""
    offset = np.abs((index - peak_time) / pd.Timedelta(days=1))
    outside = offset > half_width_days
    if not outside.any():
        return {"outside_max_change_pct": 0.0, "outside_vol_change_pct": 0.0,
                "outside_hours_changed": 0}
    o_flow = np.asarray(flow)[outside]
    o_scaled = np.asarray(scaled)[outside]
    good = o_flow > 0
    rel = np.abs(o_scaled[good] / o_flow[good] - 1.0) if good.any() else np.array([0.0])
    total = o_flow.sum()
    return {"outside_max_change_pct": 100.0 * float(rel.max()),
            "outside_vol_change_pct": (100.0 * (o_scaled.sum() - total) / total
                                       if total > 0 else 0.0),
            "outside_hours_changed": int((rel > 0.01).sum())}


def scale_volume_matched(flow, index, target_peak, target_vol5):
    """Flow-dependent multiplier hitting the peak exactly and the 5-day by iteration.

    f(Q) = f_out + (f_peak - f_out) * w,  w = (Q - Qmin) / (Qmax - Qmin)

    Monotone in Q, so the peak cannot be overtaken by a shoulder and the
    hydrograph keeps its ordering. f_out is solved from the 5-day volume over
    +/- VOLUME_HALF_WIDTH_DAYS, then refined: reshaping moves which 5-day window
    is the maximum, so the closed-form answer is only a first guess.

    Outside the volume window the multiplier is blended back to 1.0 over
    OUTSIDE_RETURN_DAYS, so the member rejoins the observed hydrograph instead
    of carrying the scaling through the whole lead-in and recession. The
    envelope is applied INSIDE the loop, so the 5-day is solved against the
    hydrograph that actually gets written.
    """
    peak_time = index[int(np.argmax(flow))]
    offset_days = (index - peak_time) / pd.Timedelta(days=1)
    inside = np.abs(offset_days) <= VOLUME_HALF_WIDTH_DAYS
    envelope = return_envelope(index, peak_time, VOLUME_HALF_WIDTH_DAYS)
    span = flow.max() - flow.min()
    weight = (flow - flow.min()) / span if span > 0 else np.ones_like(flow)
    f_peak = target_peak / flow.max()

    hi = float((flow * weight)[inside].sum())
    lo = float((flow * (1.0 - weight))[inside].sum())
    f_out = ((target_vol5 * VOLUME_HALF_WIDTH_DAYS * 2 * 24 - f_peak * hi) / lo
             if lo > 0 else f_peak)

    best = None
    for _ in range(VOLUME_ITERATIONS):
        factor = f_out + (f_peak - f_out) * weight
        scaled = apply_multiplier(flow, factor, envelope)
        got5 = rolling_max_mean(scaled, index, 120)
        err = got5 / target_vol5 - 1.0
        if best is None or abs(err) < abs(best[2]):
            best = (f_out, scaled, err)
        if abs(err) <= VOLUME_TOLERANCE:
            break
        # f_out moves the shoulders; nudge it against the residual
        step = (target_vol5 - got5) / max(got5, 1.0)
        f_out = f_out * (1.0 + 0.8 * step) if f_out > 0 else f_out + 0.05 * step
    f_out, scaled, err = best
    # peak is exact by construction, but correct after iteration -- through the
    # multiplier, so the returned-to region stays exactly observed
    factor = f_out + (f_peak - f_out) * weight
    correction = peak_correction(f_peak, target_peak / scaled.max())
    if abs(correction - 1.0) > 1e-12:
        scaled = apply_multiplier(flow, 1.0 + (factor - 1.0) * correction, envelope)
    info = {"f_peak": f_peak, "f_out": f_out,
            "shape_strain": f_out / f_peak if f_peak else np.nan,
            "vol5_err": err,
            "vol5_got": rolling_max_mean(scaled, index, 120)}
    info.update(outside_change(scaled, flow, index, peak_time,
                               VOLUME_HALF_WIDTH_DAYS))
    return scaled, info


def scale_linear_taper(flow, index, target_peak):
    """Fallback: multiplier is f_peak at the peak, falling linearly to 1.0.

    Everything beyond +/- TAPER_HALF_WIDTH_DAYS is left untouched. Meets the
    peak target only; the 5-day volume lands wherever it lands.
    """
    peak_time = index[int(np.argmax(flow))]
    offset_days = np.abs((index - peak_time) / pd.Timedelta(days=1))
    ramp = np.clip(1.0 - offset_days / TAPER_HALF_WIDTH_DAYS, 0.0, 1.0)
    f_peak = target_peak / flow.max()
    scaled = (1.0 + (f_peak - 1.0) * ramp) * flow
    info = {"f_peak": f_peak, "f_out": 1.0,
            "shape_strain": 1.0 / f_peak if f_peak else np.nan,
            "vol5_err": np.nan,
            "vol5_got": rolling_max_mean(scaled, index, 120)}
    info.update(outside_change(scaled, flow, index, peak_time,
                               VOLUME_HALF_WIDTH_DAYS))
    return scaled, info


def median_pool_by_dayofyear(dss_file, pathname, min_years):
    """Median simulated pool for each calendar day from the POR run."""
    series = read_dss_series(dss_file, pathname).dropna()
    frame = pd.DataFrame({"elev": series.values, "month": series.index.month,
                          "day": series.index.day, "year": series.index.year})
    daily = frame.groupby(["year", "month", "day"])["elev"].mean().reset_index()
    grouped = daily.groupby(["month", "day"])["elev"]
    median, counts = grouped.median(), grouped.count()
    return median[counts >= min_years]


def pool_on_date(table, when):
    key = (int(when.month), int(when.day))
    if key in table.index:
        return float(table.loc[key])
    if key == (2, 29) and (2, 28) in table.index:
        return float(table.loc[(2, 28)])
    return np.nan


def report_source_events(targets):
    """Print the enabled events with the stretch each one has to survive.

    The scale factor is the number that decides whether a source event is a
    reasonable basis for a given magnitude. An event already near the target
    barely has to be touched; one at half the target is being asked to supply a
    flood volume that storm never produced.
    """
    enabled = [row for row in SOURCE_EVENT_CATALOG if row[0]]
    if not enabled:
        raise SystemExit("SOURCE_EVENT_CATALOG has no events switched on.")
    small = [r for r in enabled if r[3] < MIN_SOURCE_PEAK_CFS]
    if small:
        raise SystemExit(
            "These events are below MIN_SOURCE_PEAK_CFS (%s cfs): %s.\n"
            "Hydrograph shape goes with size, so scaling a small storm up to a "
            "rare magnitude routes a shape that never belonged to a flood that "
            "big. Lower the setting deliberately if you want them anyway."
            % (format(int(MIN_SOURCE_PEAK_CFS), ","),
               ", ".join(r[1] for r in small)))

    # targets maps label -> {"peak", "vol5", "aep"}, so the largest PEAK has to
    # be pulled out of the nested dicts rather than taken off .values() directly.
    peaks = [row["peak"] for row in targets.values()
             if np.isfinite(row.get("peak", np.nan))]
    biggest = max(peaks) if peaks else np.nan
    print("\nSOURCE EVENTS  (%d of %d in the catalog switched on)"
          % (len(enabled), len(SOURCE_EVENT_CATALOG)))
    print("   %-9s %12s %7s %8s  %s"
          % ("event", "unreg peak", "ratio", "max x", "note"))
    stretched = []
    for _, label, date, peak, ratio, note in enabled:
        factor = biggest / peak if peak else np.nan
        flag = "  <-- stretched" if factor > SCALE_FACTOR_WARN else ""
        if factor > SCALE_FACTOR_WARN:
            stretched.append(label)
        print("   %-9s %12s %7.2f %8.2f  %s%s"
              % (label, format(peak, ","), ratio, factor, note, flag))
    ratios = [r[4] for r in enabled]
    print("   attenuation ratio spanned: %.2f to %.2f across %d events"
          % (min(ratios), max(ratios), len(enabled)))
    if max(ratios) < 0.9:
        print("   NOTE: no pass-through event enabled. The top of the regulated")
        print("   curve is set by storms that exhaust storage, so leaving those")
        print("   out will bias it low.")
    if stretched:
        print("   NOTE: %s must be scaled past %.1fx to reach the largest target."
              % (", ".join(stretched), SCALE_FACTOR_WARN))
        print("   The routing is fine arithmetically; the shape may not belong")
        print("   to a flood that size.")
    return enabled


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


def to_dss_values(values):
    filled = pd.Series(values).ffill().bfill()
    return np.where(np.isfinite(filled.values), filled.values, SENTINEL)


def fmt_dss(stamp):
    return pd.Timestamp(stamp).strftime("%d%b%Y %H%M").upper()


def d_part(start, n_hours):
    return "%s - %s" % (fmt_dss(start),
                        fmt_dss(pd.Timestamp(start) + pd.Timedelta(hours=n_hours - 1)))


def plot_events(events, mapping, stem):
    """Source hydrographs and the scaled family built from each."""
    # The grid follows the number of ENABLED source events. It was hard coded to
    # 2x2 back when the catalog was four events and switching more on crashed
    # here, after the DSS and the mapping CSV had already been written.
    n = len(events)
    ncols = 1 if n == 1 else (2 if n <= 8 else 3)
    nrows = int(np.ceil(n / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(6.75 * ncols, 4.5 * nrows))
    for k, ev in enumerate(events):
        ax = axes[k // ncols][k % ncols]
        hours = np.arange(len(ev["total"])) / 24.0
        for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
            scaled, _ = (scale_volume_matched(ev["total"].values, ev["index"],
                                              row["target_unreg_peak_cfs"],
                                              row["target_unreg_5day_cfs"])
                         if row["scaling_method"] == "volume_matched"
                         else scale_linear_taper(ev["total"].values, ev["index"],
                                                 row["target_unreg_peak_cfs"]))
            ax.plot(hours, scaled, lw=1.2,
                    label="%s  strain %.2f" % (row["target"], row["shape_strain"]))
        ax.plot(hours, ev["total"].values, color="k", lw=1.8, label="observed x1.00")
        ax.plot(hours, ev["cas"].fillna(0.0).values, color="0.6", lw=0.9, ls=":",
                label="local (observed)")
        # mark where the targets bind and where the member rejoins the observed
        peak_day = float(np.argmax(ev["total"].values)) / 24.0
        ax.axvspan(peak_day - VOLUME_HALF_WIDTH_DAYS,
                   peak_day + VOLUME_HALF_WIDTH_DAYS,
                   color="#2c7fb8", alpha=0.10, zorder=0)
        if APPLY_OUTSIDE_RETURN:
            after = (OUTSIDE_RETURN_DAYS if OUTSIDE_RETURN_DAYS_AFTER is None
                     else OUTSIDE_RETURN_DAYS_AFTER)
            for a, b in ((peak_day - VOLUME_HALF_WIDTH_DAYS - OUTSIDE_RETURN_DAYS,
                          peak_day - VOLUME_HALF_WIDTH_DAYS),
                         (peak_day + VOLUME_HALF_WIDTH_DAYS,
                          peak_day + VOLUME_HALF_WIDTH_DAYS + after)):
                ax.axvspan(a, b, color="#e67e22", alpha=0.10, zorder=0)
        ax.set_title("%s   %s" % (ev["label"], ev["note"]), fontsize=10)
        ax.set_xlabel("Days into the member window", fontsize=8)
        ax.set_ylabel("Unregulated flow (cfs)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)
    fig.suptitle("Synthetic source events and their scaled family "
                 "(Mossyrock inflow + Castle Rock local)\n"
                 "blue band = the +/- %.1f day volume window, orange = the "
                 "return to the observed hydrograph"
                 % VOLUME_HALF_WIDTH_DAYS, fontsize=11)
    for k in range(n, nrows * ncols):          # blank any unused panel
        axes[k // ncols][k % ncols].axis("off")
    fig.tight_layout(rect=[0, 0, 1, 1.0 - 0.16 / nrows])
    fig.savefig("%s_events.png" % stem, dpi=150)
    plt.close(fig)


def read_sim_inflows(sim_dss):
    """Mossyrock inflow and Castle Rock local as ResSim read them.

    One continuous record each, spanning every synthetic water year, not one
    record per member -- the members are separated by their dates. Returns
    (mos, cas) or (None, None) if the file or the records are not there.
    """
    if not os.path.isfile(sim_dss):
        print("As-run plots: %s not found -- skipped." % sim_dss)
        print("   The _events.png set shows the ensemble AS BUILT, which is "
              "not what")
        print("   was routed if the Dec1977 / Nov1986 members were chopped by "
              "hand.")
        return None, None
    try:
        mos = read_dss_series(sim_dss, "//%s/%s//1HOUR/%s/"
                              % (SIM_PART_B_MOS, SIM_PART_C_MOS, SIM_PART_F))
        cas = read_dss_series(sim_dss, "//%s/%s//1HOUR/%s/"
                              % (SIM_PART_B_CAS, SIM_PART_C_CAS, SIM_PART_F))
    except Exception as exc:
        print("As-run plots: could not read %s (%s) -- skipped."
              % (os.path.basename(sim_dss), exc))
        print("   Check SIM_PART_B_* / SIM_PART_F against the pathnames in the "
              "file.")
        return None, None
    return mos, cas


def sim_member_total(mos, cas, row, offset_hours=0):
    """One member's routed inflow total, on hours-into-the-window.

    The member is cut out of the continuous record by its own dates, the same
    way #Adjusted_Peak_Record.py cuts events out of a period-of-record run.
    Returns (hours, total) or (None, None) when the window is not in the file.

    offset_hours corrects the stamping convention between the ensemble as
    written and the simulation as written back -- see detect_sim_hour_offset.
    """
    start = pd.Timestamp(row["real_start"]) - pd.Timedelta(hours=offset_hours)
    end = pd.Timestamp(row["real_end"]) - pd.Timedelta(hours=offset_hours)
    index = pd.date_range(start, end, freq="h")
    # Sum with NaN PROPAGATING, then drop. Filling missing with zero instead
    # puts a spurious plunge to 0 on the plot: the simulation records carry
    # NaN at the member boundaries, so every panel ended with a vertical drop
    # that looked like the hydrograph collapsing rather than the record simply
    # stopping.
    total = (mos.reindex(index) + cas.reindex(index)).dropna()
    if not len(total):
        return None, None
    hours = (total.index - start) / pd.Timedelta(hours=1) / 24.0
    return np.asarray(hours, dtype=float), total.values


def residual_against_built(hours, total, built):
    """Compare a simulated member to its built curve, aligned by TIME.

    Not by array position. sim_member_total drops the missing hours at the
    member boundaries, so position i in the simulated array is not hour i --
    comparing positionally slides the whole record and reports differences
    everywhere, which is how the first version of this made all 48 members look
    hand-edited.

    Returns (max_abs, index_into_built_of_max, run_value, built_value) or None.
    """
    if hours is None or total is None or not len(total):
        return None
    hour_index = np.rint(np.asarray(hours) * 24.0).astype(int)
    inside = (hour_index >= 0) & (hour_index < len(built))
    if inside.sum() < 100:
        return None
    hi = hour_index[inside]
    run = np.asarray(total)[inside]
    ref = np.asarray(built)[hi]
    diff = run - ref
    j = int(np.nanargmax(np.abs(diff)))
    return float(abs(diff[j])), int(hi[j]), float(run[j]), float(ref[j])


def detect_sim_hour_offset(events, mapping, mos, cas):
    """Hours to shift the simulation record to line up with the build.

    DSS stamps are END of period, and the ensemble is written one step after
    the label start while read_dss_series shifts back to period-beginning. The
    two conventions do not cancel, and the simulation comes back one hour out.
    Uncorrected, every member is compared against its own curve shifted by an
    hour, which produces large apparent differences on the steep limbs while
    the peaks still match -- easily mistaken for a hand edit on all 48 members.

    Detected rather than assumed: each member is tested over a range of shifts
    and the one that minimises the largest residual is kept, then the MODE
    across members is adopted. The mode is what makes this safe -- the members
    that were edited by hand disagree with their built curve at every shift, so
    letting them vote individually would be wrong, but they are outnumbered.
    """
    if SIM_HOUR_OFFSET is not None:
        print("As-run plots: hour offset forced to %+d" % SIM_HOUR_OFFSET)
        return SIM_HOUR_OFFSET
    votes = []
    for ev in events:
        for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
            built, _ = scale_volume_matched(
                ev["total"].values, ev["index"],
                row["target_unreg_peak_cfs"], row["target_unreg_5day_cfs"])
            best, best_resid = 0, None
            for shift in range(-SIM_OFFSET_SEARCH, SIM_OFFSET_SEARCH + 1):
                hours, total = sim_member_total(mos, cas, row, shift)
                if total is None:
                    continue
                resid = residual_against_built(hours, total, built)
                if resid is None:
                    continue
                if best_resid is None or resid[0] < best_resid:
                    best, best_resid = shift, resid[0]
            if best_resid is not None:
                votes.append(best)
    if not votes:
        return 0
    offset = int(pd.Series(votes).mode().iloc[0])
    agree = votes.count(offset)
    print("As-run plots: simulation is %+d hour(s) off the build index "
          "(%d of %d members agree)" % (offset, agree, len(votes)))
    if agree < 0.6 * len(votes):
        print("   *** less than 60%% agreement -- the offset is not clean. "
              "Check the")
        print("   mapping CSV against the simulation before trusting these "
              "plots.")
    return offset


def plot_events_overlay(events, stem):
    """Every source storm on one pair of axes, unscaled.

    This is the figure that answers "did you look at a range of shapes?" --
    the observed hydrographs, no scaling, aligned on their own peak.
    """
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for ev in events:
        v = ev["total"].values
        hours = (np.arange(len(v)) - float(np.argmax(v))) / 24.0
        ax.plot(hours, v, lw=1.4, label="%s" % ev["label"])
    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel("Days from the observed peak")
    ax.set_ylabel("Unregulated flow (cfs)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: format(int(v), ",")))
    ax.set_title("The %d observed source storms, unscaled\n"
                 "Mossyrock inflow + Castle Rock local, aligned on each "
                 "storm's own peak" % len(events), fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig("%s_events_overlay.png" % stem, dpi=150)
    print("Wrote %s_events_overlay.png" % stem)


def plot_event_single(events, mapping, stem, label=None):
    """One source storm and the whole family scaled from it."""
    pick = None
    if label:
        pick = next((e for e in events if e["label"] == label), None)
    if pick is None:
        pick = events[0]
    ev = pick
    fig, ax = plt.subplots(figsize=(10, 6.5))
    hours = np.arange(len(ev["total"])) / 24.0
    for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
        scaled, _ = (scale_volume_matched(ev["total"].values, ev["index"],
                                          row["target_unreg_peak_cfs"],
                                          row["target_unreg_5day_cfs"])
                     if row["scaling_method"] == "volume_matched"
                     else scale_linear_taper(ev["total"].values, ev["index"],
                                             row["target_unreg_peak_cfs"]))
        ax.plot(hours, scaled, lw=1.6,
                label="%s  (strain %.2f)" % (row["target"], row["shape_strain"]))
    ax.plot(hours, ev["total"].values, color="k", lw=2.2,
            label="observed x1.00")
    peak_day = float(np.argmax(ev["total"].values)) / 24.0
    ax.axvspan(peak_day - VOLUME_HALF_WIDTH_DAYS,
               peak_day + VOLUME_HALF_WIDTH_DAYS,
               color="#2c7fb8", alpha=0.12, zorder=0,
               label="+/- %.1f day volume window" % VOLUME_HALF_WIDTH_DAYS)
    ax.set_xlabel("Days into the member window")
    ax.set_ylabel("Unregulated flow (cfs)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: format(int(v), ",")))
    ax.set_title("Scaling one source storm: %s   %s\n"
                 "each member is scaled to its own target peak and 5-day "
                 "volume" % (ev["label"], ev["note"]), fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("%s_event_single.png" % stem, dpi=150)
    print("Wrote %s_event_single.png  (source event %s)" % (stem, ev["label"]))


def plot_events_asrun(events, mapping, stem, mos, cas, offset_hours=0):
    """The same panels as plot_events, drawn from the simulation's own inflows.

    This is the honest picture of what was routed: it carries the hand-edited
    Dec1977 and Nov1986 members rather than the two-peaked versions the scaling
    functions produce.
    """
    n = len(events)
    ncols = 1 if n == 1 else (2 if n <= 8 else 3)
    nrows = int(np.ceil(n / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(6.75 * ncols, 4.5 * nrows))
    drawn, missing = 0, []
    for k, ev in enumerate(events):
        ax = axes[k // ncols][k % ncols]
        obs_hours = np.arange(len(ev["total"])) / 24.0
        for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
            hours, total = sim_member_total(mos, cas, row, offset_hours)
            if hours is None:
                missing.append("%s %s" % (row["event"], row["target"]))
                continue
            line = ax.plot(hours, total, lw=1.2,
                           label="%s  peak %s" % (row["target"],
                                                  format(int(total.max()), ",")))[0]
            drawn += 1
            if ASRUN_SHOW_BUILT:
                built, _ = scale_volume_matched(
                    ev["total"].values, ev["index"],
                    row["target_unreg_peak_cfs"], row["target_unreg_5day_cfs"])
                ax.plot(obs_hours, built, lw=0.8, ls="--",
                        color=line.get_color(), alpha=0.8)
        ax.plot(obs_hours, ev["total"].values, color="k", lw=1.8,
                label="observed x1.00")
        ax.plot(obs_hours, ev["cas"].fillna(0.0).values, color="0.6", lw=0.9,
                ls=":", label="local (observed)")
        peak_day = float(np.argmax(ev["total"].values)) / 24.0
        ax.axvspan(peak_day - VOLUME_HALF_WIDTH_DAYS,
                   peak_day + VOLUME_HALF_WIDTH_DAYS,
                   color="#2c7fb8", alpha=0.10, zorder=0)
        ax.set_title("%s   %s" % (ev["label"], ev["note"]), fontsize=10)
        ax.set_xlabel("Days into the member window", fontsize=8)
        ax.set_ylabel("Unregulated inflow (cfs)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)
    fig.suptitle("Synthetic inflows AS RUN -- read back from %s\n"
                 "Mossyrock inflow + Castle Rock local, including any members "
                 "edited by hand after the ensemble was written%s"
                 % (os.path.basename(SIM_DSS),
                    "; dashed = as built" if ASRUN_SHOW_BUILT else ""),
                 fontsize=11)
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    fig.tight_layout(rect=[0, 0, 1, 1.0 - 0.16 / nrows])
    fig.savefig("%s_events_asrun.png" % stem, dpi=150)
    plt.close(fig)
    print("As-run plots: %d member(s) drawn from %s"
          % (drawn, os.path.basename(SIM_DSS)))
    if missing:
        print("   %d member(s) had no window in the simulation: %s"
              % (len(missing), ", ".join(missing[:8])))
        print("   Those members were not run, or the mapping CSV does not "
              "match the simulation.")
    return drawn


def compare_built_to_asrun(events, mapping, mos, cas, offset_hours=0):
    """Report which members differ between the built ensemble and the run.

    This is what makes a hand edit visible instead of inferred. A member that
    was chopped shows up as a lower pre-peak maximum against an unchanged peak;
    a member nobody touched matches to within rounding.
    """
    rows = []
    for ev in events:
        for _, row in mapping[mapping["event"] == ev["label"]].iterrows():
            hours, total = sim_member_total(mos, cas, row, offset_hours)
            if hours is None:
                continue
            built, _ = scale_volume_matched(
                ev["total"].values, ev["index"],
                row["target_unreg_peak_cfs"], row["target_unreg_5day_cfs"])
            # Largest hourly difference anywhere in the window, aligned by
            # time. Looking only at the pre-peak MAXIMUM missed Dec1977
            # entirely: its chop lowered a shoulder without moving that
            # shoulder's maximum enough to trip a threshold.
            found = residual_against_built(hours, total, built)
            if found is None:
                continue
            max_abs, hour_at, run_there, built_there = found
            rows.append({"event": row["event"], "target": row["target"],
                         "peak_built": float(np.max(built)),
                         "peak_run": float(np.max(total)),
                         "at_day": hour_at / 24.0,
                         "built_there": built_there,
                         "run_there": run_there,
                         "max_abs_diff": max_abs})
    if not rows:
        return
    table = pd.DataFrame(rows)
    table["diff_pct_of_peak"] = 100.0 * table["max_abs_diff"] / table["peak_built"]
    edited = table[table["diff_pct_of_peak"] > EDIT_DETECT_PCT]
    print("\nBUILT vs AS RUN  (largest hourly difference, as %% of the member "
          "peak)")
    print("   peaks match to within %.1f%% on every member"
          % (100 * (table["peak_run"] / table["peak_built"] - 1.0).abs().max()))
    if not len(edited):
        print("   nothing differs by more than %.0f%% -- no hand edit detected"
              % EDIT_DETECT_PCT)
        return
    print("   %d member(s) differ by more than %.0f%%:"
          % (len(edited), EDIT_DETECT_PCT))
    for _, r in edited.sort_values("diff_pct_of_peak", ascending=False).iterrows():
        print("      %-8s %-7s day %5.2f  %8s -> %8s  (%+.0f cfs, %.1f%% of peak)"
              % (r["event"], r["target"], r["at_day"],
                 format(int(r["built_there"]), ","),
                 format(int(r["run_there"]), ","),
                 r["run_there"] - r["built_there"], r["diff_pct_of_peak"]))
    print("   That is the hand edit, and it is why the as-run set is the one "
          "to publish.")


def plot_pools(events, bases, highest_pick, stem):
    """Starting pool for every ENABLED basis, per event date."""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(events))
    width = min(0.26, 0.8 / max(len(bases), 1))
    colors = {"rulecurve": "#8e44ad", "duration50": "#2c7fb8",
              "observed": "#c0392b", "median_por": "#16a085",
              "highest": "#e67e22"}
    offset0 = (len(bases) - 1) / 2.0
    for i, basis in enumerate(bases):
        vals = [ev["pools"].get(basis, np.nan) for ev in events]
        ax.bar(x + (i - offset0) * width, vals, width,
               color=colors.get(basis, "0.5"), label=basis)
    ax.set_xticks(x)
    ax.set_xticklabels([ev["label"] for ev in events])
    ax.set_ylim(700, 790)
    ax.set_ylabel("Starting pool elevation (ft)")
    off = [b for b in POOL_BASIS_ORDER if b not in bases]
    ax.set_title("Starting pool basis by source event -- each becomes its own "
                 "ensemble member\non: %s%s"
                 % (", ".join(bases),
                    "     off: %s" % ", ".join(off) if off else ""))
    if "highest" in bases and highest_pick:
        ax.text(0.995, 0.02,
                "highest picked: %s" % ", ".join("%s=%s" % (k, v) for k, v
                                                 in highest_pick.items()),
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color=colors["highest"])
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig("%s_pools.png" % stem, dpi=150)
    plt.close(fig)


def main():
    for path in (os.path.dirname(MAPPING_CSV), os.path.dirname(PLOT_STEM)):
        if path and not os.path.isdir(path):
            os.makedirs(path)

    mos = read_dss_series(IN_DSS, PATH_MOS_IN).dropna()
    cas = read_dss_series(IN_DSS, PATH_CAS_LOCAL).dropna()
    if CLIP_NEGATIVE_FLOW:
        mos = mos.clip(lower=0.0)
        cas = cas.clip(lower=0.0)
    elev_daily = read_dss_series(OBS_DSS, PATH_MOS_ELEV_DAILY).dropna()
    elev_hourly = daily_to_hourly(elev_daily, ELEV_DAILY_TO_HOURLY)
    dur50 = duration_50_by_dayofyear(elev_daily, DURATION_50_MIN_YEARS)

    n_hours = (WINDOW_BEFORE_DAYS + WINDOW_AFTER_DAYS) * 24
    elev_hours = n_hours + LOOKBACK_DAYS * 24
    ens_start = pd.Timestamp(ENS_LABEL_START)
    elev_ens_start = ens_start - pd.Timedelta(days=LOOKBACK_DAYS)
    flow_d = d_part(ens_start, n_hours)
    elev_d = d_part(elev_ens_start, elev_hours)
    flow_write = (ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()
    elev_write = (elev_ens_start + pd.Timedelta(hours=1)).strftime("%d%b%Y %H:%M:%S").upper()

    bases = [b for b in POOL_BASIS_ORDER if POOL_BASES_ENABLED.get(b, False)]
    if not bases:
        raise SystemExit("POOL_BASES_ENABLED has nothing switched on -- "
                         "there would be no members to write.")
    if bases == ["highest"]:
        raise SystemExit("'highest' picks among the OTHER enabled bases -- "
                         "switch at least one of those on as well.")
    if POOL_FALLBACK_BASIS not in ("rulecurve",) and POOL_FALLBACK_BASIS not in bases:
        raise SystemExit("POOL_FALLBACK_BASIS is '%s', which is not switched on "
                         "in POOL_BASES_ENABLED." % POOL_FALLBACK_BASIS)
    highest_pick = {}
    targets = read_targets(UNREG_FREQ_CSV, FREQ_VALUE_COL)
    por_medians = None
    if "median_por" in bases:
        if os.path.isfile(POR_ELEV_DSS):
            por_medians = median_pool_by_dayofyear(POR_ELEV_DSS, POR_ELEV_PATH,
                                                   POR_ELEV_MIN_YEARS)
            print("Median POR pool: %d calendar days from %s"
                  % (len(por_medians), POR_ELEV_DSS))
        else:
            print("'median_por' is switched on but the POR run was not found:")
            print("   %s" % POR_ELEV_DSS)
            print("   Falling back to %s; every member is tagged so."
                  % POOL_FALLBACK_BASIS)

    print("=" * 78)
    print("Window     : %d days before / %d days after the peak (%d hours)"
          % (WINDOW_BEFORE_DAYS, WINDOW_AFTER_DAYS, n_hours))
    print("Design     : %d events x %d magnitudes x %d pool bases = %d members"
          % (len(SOURCE_EVENTS), len(TARGET_AEPS) + 1, len(bases),
             len(SOURCE_EVENTS) * (len(TARGET_AEPS) + 1) * len(bases)))
    print("Scaling    : %s -- both records scaled, hourly proportion preserved"
          % SCALING_METHOD)
    print("Start pool : %d of %d bases ON"
          % (len(bases), len(POOL_BASIS_ORDER)))
    for basis in POOL_BASIS_ORDER:
        on = POOL_BASES_ENABLED.get(basis, False)
        print("             [%s] %-11s %s"
              % ("x" if on else " ", basis,
                 "" if on else "(off -- no members written)"))
    print("Output     : %s" % OUT_DSS)
    print("=" * 78)
    report_source_events(targets)

    events = []
    for label, peak_date, note in SOURCE_EVENTS:
        centre = pd.Timestamp(peak_date)
        start = (centre - pd.Timedelta(days=WINDOW_BEFORE_DAYS)).normalize()
        index = pd.date_range(start, periods=n_hours, freq="h")
        mos_w = mos.reindex(index)
        cas_w = cas.reindex(index)
        total = mos_w.fillna(0.0) + cas_w.fillna(0.0)
        elev_index = pd.date_range(start - pd.Timedelta(days=LOOKBACK_DAYS),
                                   periods=elev_hours, freq="h")
        obs_pool = elev_hourly.reindex(elev_index)
        rc_pool = pd.Series(rule_curve_on_index(elev_index), index=elev_index)
        pools = {"rulecurve": float(rc_pool.loc[start])}
        series_by_basis = {"rulecurve": rc_pool}

        if "duration50" in bases:
            series, value = pool_series_for(dur50, elev_index, start,
                                            POOL_SERIES_STYLE)
            pools["duration50"] = value
            series_by_basis["duration50"] = series if np.isfinite(value) else rc_pool

        if "observed" in bases:
            value = (float(obs_pool.loc[start])
                     if start in obs_pool.index and np.isfinite(obs_pool.loc[start])
                     else np.nan)
            pools["observed"] = value
            # full-length lookback record, same convention as #Create_ObsRC_Ensembles
            series_by_basis["observed"] = obs_pool if np.isfinite(value) else rc_pool

        if "median_por" in bases:
            if por_medians is not None:
                series, value = pool_series_for(por_medians, elev_index, start,
                                                POOL_SERIES_STYLE)
            else:
                series, value = rc_pool, np.nan
            pools["median_por"] = value
            series_by_basis["median_por"] = series if np.isfinite(value) else rc_pool

        if "highest" in bases:
            # the highest of the OTHER enabled bases that actually have a value
            candidates = {k: pools[k] for k in bases
                          if k != "highest" and np.isfinite(pools.get(k, np.nan))}
            if candidates:
                pick = max(candidates, key=candidates.get)
                pools["highest"] = candidates[pick]
                series_by_basis["highest"] = series_by_basis[pick]
                highest_pick[label] = pick
            else:
                pools["highest"] = np.nan
                series_by_basis["highest"] = rc_pool
                highest_pick[label] = "none"
        events.append({"label": label, "note": note, "peak_date": centre,
                       "start": start, "index": index, "elev_index": elev_index,
                       "mos": mos_w, "cas": cas_w, "total": total,
                       "obs_peak": float(total.max()),
                       "local_share": float(cas_w.loc[total.idxmax()] /
                                            total.max()) if total.max() > 0 else np.nan,
                       "pools": pools, "pool_series": series_by_basis,
                       "missing": int(mos_w.isna().sum() + cas_w.isna().sum())})

    def build_set(out_dss, method):
      mapping_rows = []
      member = 0
      # HecDss.Open() does not truncate an existing file -- it opens the file
      # in place and only overwrites the pathnames this run actually writes.
      # Every design change here (event catalog, pool bases, member count) has
      # moved which member number gets which event, so a stale out_dss left
      # over from an earlier run keeps old members' data under F-parts this
      # run no longer writes, and CAN leave stale storage blocks behind a
      # pathname this run DOES write if the new record is shorter than the
      # old one at that pathname. Either way the DSS and the mapping CSV it is
      # written alongside stop agreeing about what a given member number is.
      # Delete first so every run starts from an empty file, same as
      # #Adjusted_Peak_Record.py does for its own output DSS.
      #
      # THIS IS ALSO WHAT DESTROYS THE HAND-CHOPPED Dec1977/Nov1986 MEMBERS.
      # See the warning at the top of this file: the mini peaks are removed by
      # hand in DSSVue after this script runs, and deleting the file throws
      # that away silently. The chop has to be redone before the ResSim run.
      if os.path.exists(out_dss):
          print("   NOTE: %s existed and is being replaced. Any hand edits in it"
                % os.path.basename(out_dss))
          print("         (the Dec1977 / Nov1986 mini-peak chops) are now gone")
          print("         and must be redone before the ResSim run.")
          os.remove(out_dss)
      with HecDss.Open(out_dss, version=OUT_DSS_VERSION) as dst:
          for ev in events:
              for target_label, target in targets.items():
                  target_peak, target_vol5 = target["peak"], target["vol5"]
                  if method == "volume_matched":
                      tot_scaled, info = scale_volume_matched(
                          ev["total"].values, ev["index"], target_peak, target_vol5)
                  else:
                      tot_scaled, info = scale_linear_taper(
                          ev["total"].values, ev["index"], target_peak)
                  # split the scaled total back onto the two records in the
                  # observed proportion at each hour, so the coincidence between
                  # reservoir inflow and local is preserved exactly
                  share = np.divide(ev["mos"].fillna(0.0).values,
                                    np.maximum(ev["total"].values, 1e-9))
                  mos_scaled = tot_scaled * share
                  cas_scaled = tot_scaled * (1.0 - share)
                  factor = float(np.nanmedian(np.divide(
                      tot_scaled, np.maximum(ev["total"].values, 1e-9))))
                  for basis in bases:
                      pool = ev["pools"].get(basis, np.nan)
                      basis_used = basis
                      if not np.isfinite(pool):
                          pool = ev["pools"][POOL_FALLBACK_BASIS]
                          basis_used = "%s->%s (unavailable)" % (
                              basis, POOL_FALLBACK_BASIS)
                      member += 1
                      f_part = "C:%06d|%s" % (member, ENS_SUFFIX)
                      mos_v = to_dss_values(mos_scaled)
                      cas_v = to_dss_values(cas_scaled)
                      pool_series = ev["pool_series"][basis]
                      if basis_used != basis:
                          pool_series = ev["pool_series"][POOL_FALLBACK_BASIS]
                      elev_v = to_dss_values(pool_series.values)
                      rc_v = to_dss_values(rule_curve_on_index(ev["elev_index"]))
                      for parts, vals, units, dpart, wstart in [
                              (("", "MOSSYROCK", "FLOW-IN"), mos_v, "CFS", flow_d, flow_write),
                              (("", "CASTLE ROCK", "FLOW-LOCAL"), cas_v, "CFS", flow_d, flow_write),
                              (("", "MOS", "ELEV"), elev_v, "FEET", elev_d, elev_write),
                              (("", "MOS", "ELEV-RULECURVE"), rc_v, "FEET", elev_d, elev_write)]:
                          pathname = "/%s/%s/%s/%s/1HOUR/%s/" % (
                              parts[0], parts[1], parts[2], dpart, f_part)
                          dst.put_ts(build_container(pathname, vals, wstart, units,
                                                     "INST-VAL", 60))
                      # Shift by WATER year, not calendar year. Shifting the
                      # calendar year puts an Oct-Dec event in synthetic WY
                      # synth_year+1 while a Jan-Sep event lands in synth_year,
                      # so autumn and winter members collide in one water year
                      # and leave others empty (and 30-day windows that cross
                      # 01 Jan can overlap the next member outright).
                      synth_year = SYNTH_YEAR_BASE + member
                      synth_start = ev["start"] + pd.DateOffset(
                          years=synth_year - water_year(ev["start"]))
                      mapping_rows.append({
                          "member": member, "ensemble_f_part": f_part,
                          "event": ev["label"], "event_note": ev["note"],
                          "target": target_label, "target_aep": target["aep"],
                          "scaling_method": method,
                          "target_unreg_peak_cfs": round(target_peak, 1),
                          "target_unreg_5day_cfs": round(target_vol5, 1),
                          "f_peak": round(info["f_peak"], 4),
                          "f_out": round(info["f_out"], 4),
                          "shape_strain": round(info["shape_strain"], 3),
                          "scaled_unreg_peak_cfs": round(float(tot_scaled.max()), 1),
                          "scaled_unreg_5day_cfs": round(info["vol5_got"], 1),
                          "vol5_error_pct": (round(100 * info["vol5_err"], 2)
                                             if np.isfinite(info["vol5_err"]) else np.nan),
                          "outside_max_change_pct": round(
                              info["outside_max_change_pct"], 2),
                          "outside_vol_change_pct": round(
                              info["outside_vol_change_pct"], 2),
                          "outside_hours_changed": info["outside_hours_changed"],
                          "scale_factor": round(factor, 4),
                          "observed_unreg_peak_cfs": round(ev["obs_peak"], 1),
                          "local_share_at_peak": round(ev["local_share"], 4),
                          "pool_basis": basis,
                          "pool_basis_used": basis_used,
                          "start_pool_ft": round(pool, 2),
                          "synth_water_year": synth_year,
                          "real_start": synth_start,
                          "real_end": synth_start + pd.Timedelta(hours=n_hours - 1),
                          "source_start": ev["start"],
                          "source_peak_date": ev["peak_date"],
                          "ensemble_start": ens_start, "hours": n_hours,
                          "elev_ensemble_start": elev_ens_start,
                          "elev_hours": elev_hours,
                          "missing_input_hours": ev["missing"]})

      return pd.DataFrame(mapping_rows)

    if PLOTS_ONLY:
        # No DSS is written or deleted in this mode -- see PLOTS_ONLY above.
        if not os.path.isfile(MAPPING_CSV):
            raise SystemExit(
                "PLOTS_ONLY needs an existing mapping CSV and there is none "
                "at\n   %s\nRun once with PLOTS_ONLY = False to build the "
                "ensemble first." % MAPPING_CSV)
        mapping = pd.read_csv(MAPPING_CSV)
        print("\nPLOTS_ONLY: %d members read from %s"
              % (len(mapping), MAPPING_CSV))
        print("   The ensemble DSS was NOT touched, so any hand-edited members "
              "in it are intact.")
    else:
        mapping = build_set(OUT_DSS, SCALING_METHOD)
        mapping.to_csv(MAPPING_CSV, index=False)
        if ALSO_WRITE_LINEAR_TAPER and SCALING_METHOD != "linear_taper":
            alt = build_set(OUT_DSS_ALT, "linear_taper")
            alt.to_csv(MAPPING_CSV.replace(".csv", "_lineartaper.csv"), index=False)
            print("\nFallback set written: %s" % OUT_DSS_ALT)
            print("   peak targets met; 5-day lands at %.0f-%.0f%% of target"
                  % (100 * (alt["scaled_unreg_5day_cfs"] / alt["target_unreg_5day_cfs"]).min(),
                     100 * (alt["scaled_unreg_5day_cfs"] / alt["target_unreg_5day_cfs"]).max()))
    plot_events(events, mapping, PLOT_STEM)
    if PLOT_SUMMARY_FIGURES:
        plot_events_overlay(events, PLOT_STEM)
        plot_event_single(events, mapping, PLOT_STEM, FEATURED_EVENT)
    plot_pools(events, bases, highest_pick, PLOT_STEM)

    # --- the same panels, drawn from what ResSim actually read ---------------
    if PLOT_AS_RUN:
        sim_mos, sim_cas = read_sim_inflows(SIM_DSS)
        if sim_mos is not None:
            offset = detect_sim_hour_offset(events, mapping, sim_mos, sim_cas)
            if plot_events_asrun(events, mapping, PLOT_STEM, sim_mos, sim_cas,
                                 offset):
                compare_built_to_asrun(events, mapping, sim_mos, sim_cas,
                                       offset)

    print("\nSOURCE EVENTS")
    for ev in events:
        pool_txt = "  ".join(
            "%s %.1f" % (k, ev["pools"][k]) if np.isfinite(ev["pools"].get(k, np.nan))
            else "%s n/a" % k for k in bases)
        print("   %-8s peak %7.0f  local %4.1f%%  gaps %d  pools: %s"
              % (ev["label"], ev["obs_peak"], 100 * ev["local_share"],
                 ev["missing"], pool_txt))
    print("\nSCALING -- f_peak at the peak, f_out at low flow, and the shape")
    print("strain between them (1.0 = uniform scaling, far from 1 = reshaped)")
    show = mapping[["event", "target", "f_peak", "f_out", "shape_strain",
                    "scaled_unreg_peak_cfs", "target_unreg_peak_cfs",
                    "scaled_unreg_5day_cfs", "target_unreg_5day_cfs",
                    "vol5_error_pct"]].copy()
    print(show.round(3).to_string(index=False))
    strained = mapping[(mapping["shape_strain"] > SHAPE_STRAIN_WARN) |
                       (mapping["shape_strain"] < 1.0 / SHAPE_STRAIN_WARN)]
    if len(strained):
        print("\n   %d members reshape by more than %.1fx between peak and"
              % (len(strained), SHAPE_STRAIN_WARN))
        print("   shoulder. Sharp events pushed to a volume target get their")
        print("   shoulders inflated; sustained events pushed to a peak target")
        print("   get them cut. Inspect those hydrographs before using them:")
        for _, r in strained.iterrows():
            print("      %-8s %-7s strain %.2f" % (r["event"], r["target"],
                                                   r["shape_strain"]))
    worst = mapping["vol5_error_pct"].abs().max()
    if np.isfinite(worst):
        print("\n   5-day volume: worst miss %.2f%% (tolerance %.1f%%)"
              % (worst, 100 * VOLUME_TOLERANCE))

    print("\nOUTSIDE THE +/- %.1f DAY WINDOW" % VOLUME_HALF_WIDTH_DAYS)
    if APPLY_OUTSIDE_RETURN:
        after = (OUTSIDE_RETURN_DAYS if OUTSIDE_RETURN_DAYS_AFTER is None
                 else OUTSIDE_RETURN_DAYS_AFTER)
        print("   returning to the observed hydrograph over %.2f day(s) before "
              "and %.2f after," % (OUTSIDE_RETURN_DAYS, after))
        print("   %s shaped, so a member is back on the observed flow %.2f days "
              "past the peak." % (OUTSIDE_RETURN_SHAPE,
                                  VOLUME_HALF_WIDTH_DAYS + after))
    else:
        print("   APPLY_OUTSIDE_RETURN is False -- the multiplier runs across "
              "the whole member window.")
    sub = mapping[["event", "target", "outside_max_change_pct",
                   "outside_vol_change_pct", "outside_hours_changed"]]
    sub = sub.drop_duplicates(subset=["event", "target"])
    print("   largest single-hour change outside the window: %.2f%%"
          % sub["outside_max_change_pct"].max())
    print("   volume outside the window vs observed: %+.2f%% to %+.2f%%"
          % (sub["outside_vol_change_pct"].min(),
             sub["outside_vol_change_pct"].max()))
    print("   (hours changed by more than 1%%: max %d of %d in a member)"
          % (sub["outside_hours_changed"].max(),
             (WINDOW_BEFORE_DAYS + WINDOW_AFTER_DAYS) * 24))
    print("   Turn OUTSIDE_RETURN_DAYS down to rejoin sooner, up for a gentler "
          "join.")
    print("\nSynthetic water years %d-%d, one per member. source_start holds the"
          % (mapping["synth_water_year"].min(), mapping["synth_water_year"].max()))
    print("true event date.")
    check_synthetic_calendar(mapping)
    print("\nMembers written : %d   records: %d (4 per member)"
          % (len(mapping), len(mapping) * 4))
    print("Mapping CSV     : %s" % MAPPING_CSV)
    print("Plots           : %s_events.png (as built), %s_pools.png"
          % (PLOT_STEM, PLOT_STEM))
    if PLOT_AS_RUN and os.path.isfile(SIM_DSS):
        print("                  %s_events_asrun.png (as run -- use this one)"
              % PLOT_STEM)


main()
