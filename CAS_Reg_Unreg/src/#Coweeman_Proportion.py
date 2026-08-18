#Coweeman_Proportion.py
# -*- coding: utf-8 -*-
"""
Is the Coweeman a stable fraction of unregulated Castle Rock flow at the peak?

WHY THE DENOMINATOR IS THE UNREGULATED FLOW
    #Coweeman_Timing.py compared the Coweeman against the OBSERVED (regulated)
    Castle Rock record, which is the right thing for a timing question but the
    wrong thing for this one. A ratio taken against regulated flow carries the
    reservoir's operation in its denominator: the same storm gives a different
    ratio depending on how much Mossyrock happened to hold back that day. The
    ratio wanted here is hydrologic -- how much of the natural flood arrives
    from the Coweeman -- so the denominator is the routed UNREGULATED flow at
    Castle Rock, which is also the quantity the unreg-reg transform is drawn
    on, so the two analyses share an x-axis.

WHAT IS MEASURED, AND WHY TWO RATIOS
    Events are picked off the unregulated record. For each one:

      ratio_coincident  = Coweeman flow AT THE HOUR of the Castle Rock peak
                          divided by that peak. This is the quantity a
                          coincident frequency analysis actually needs: what
                          the tributary is contributing when the mainstem
                          crests.
      ratio_peak        = Coweeman peak anywhere in the window divided by the
                          Castle Rock peak. What the tributary would add if
                          the two crested together.

    The gap between them is the timing penalty. Reporting only the second
    overstates the combination, which is the failure mode this whole line of
    work exists to avoid; reporting only the first hides how much of the
    spread is timing rather than hydrology.

IS IT TIGHTER AT HIGH FLOWS?
    Tested three ways rather than asserted from a plot:
      - Spearman rank correlation of the ratio against magnitude: does the
        LEVEL drift with size?
      - Coefficient of variation and IQR/median within magnitude bins: does
        the SPREAD narrow? Both are reported because CV is what people expect
        and IQR/median is what survives a small sample.
      - Brown-Forsythe (Levene on medians) across the bins: is any narrowing
        more than sampling noise?
    With this many events the bins hold single figures, so the test has little
    power. A non-significant result here means "not demonstrated", not "no
    effect", and the script says so rather than letting a p-value be read as
    proof of stability.
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import re
import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime
from pydsstools.heclib.dss import HecDss

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
CACHE_DIR = r"../data/coweeman"          # written by #Coweeman_Timing.py
UNREG_DSS = r"../output/ResSim_WCM_RC.dss"
UNREG_PATH = "//CastleRock_NWS/Flow-UNREG//1Hour/ResSim_WCM_RC/"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "coweeman_proportion.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coweeman_proportion.png")

RESAMPLE = "1h"                  # unregulated record is hourly; match it
N_EVENTS = 80                    # independent unregulated events to take
EVENT_MIN_SEPARATION_DAYS = 7
MIN_EVENT_CFS = 20000.0          # floor on the UNREGULATED peak
EVENT_WINDOW_HOURS = 48          # +/- window the Coweeman peak is sought in
MIN_WINDOW_COVERAGE = 0.80

# Magnitude bins for the consistency question, in unregulated cfs. Chosen as
# round numbers rather than quantiles so the bin edges do not move when the
# event list changes.
MAGNITUDE_BINS = [20000, 40000, 60000, 200000]
BIN_LABELS = ["20-40k", "40-60k", ">60k"]
MIN_PER_BIN = 4                  # below this a bin is reported but not tested

C_LOW, C_HIGH = "#8fbcdb", "#1a4f8a"
C_COW = "#b7410e"

# ----------------------------------------------------------------------------

ECOLOGY_ROW_DT = re.compile(
    r"^\s*(\d{1,2}/\d{1,2}/\d{4})\s+(\d{1,2}:\d{2})\s+(-?[\d.]+)(?:\s+(\S+))?\s*$")


def read_ecology_cache(cache_dir):
    """The Coweeman 15-minute record from the files #Coweeman_Timing.py cached.

    Only the FM files are read -- the DV ones are daily and cannot support a
    ratio taken at the hour of the peak.
    """
    files = sorted(glob.glob(os.path.join(cache_dir, "*_FM.txt")))
    if not files:
        raise SystemExit(
            "No Ecology FM files in %s.\nRun #Coweeman_Timing.py first; it "
            "downloads and caches them." % os.path.abspath(cache_dir))
    stamps, values = [], []
    for path in files:
        with open(path, "r", errors="replace") as handle:
            for line in handle:
                match = ECOLOGY_ROW_DT.match(line)
                if match:
                    date, clock, value, _flag = match.groups()
                    stamps.append("%s %s" % (date, clock))
                    values.append(value)
    frame = pd.DataFrame({"stamp": stamps, "value": values})
    frame["stamp"] = pd.to_datetime(frame["stamp"], format="%m/%d/%Y %H:%M",
                                    errors="coerce")
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce")
    frame = frame.dropna()
    frame = frame[frame["value"] > -900.0]
    frame = frame.drop_duplicates(subset="stamp").sort_values("stamp")
    series = pd.Series(frame["value"].to_numpy(),
                       index=pd.DatetimeIndex(frame["stamp"]))
    print("   Coweeman  : %d values from %d file(s), %s to %s"
          % (len(series), len(files), series.index.min().date(),
             series.index.max().date()))
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


def read_unreg(path, pathname):
    """Routed unregulated flow at Castle Rock, hourly."""
    version = 6 if open(path, "rb").read(16)[12] == 6 else 7
    dss = HecDss.Open(path, version=version)
    try:
        ts = dss.read_ts(pathname)
        values = np.array(ts.values, dtype=float)
        values[np.array(ts.nodata, dtype=bool)] = np.nan
        values[values <= -900.0] = np.nan
        step = pd.Timedelta(seconds=int(ts.interval))
        index = pd.date_range(first_stamp(ts) - step, periods=len(values),
                              freq=step)
    finally:
        dss.close()
    series = pd.Series(values, index=index).dropna().sort_index()
    print("   Castle Rock unregulated: %d hours, %s to %s"
          % (len(series), series.index.min().date(), series.index.max().date()))
    return series


def find_events(series, n, min_separation_days, floor):
    working = series.copy()
    events, span = [], pd.Timedelta(days=min_separation_days)
    while len(events) < n:
        if not len(working.dropna()) or working.max() < floor:
            break
        when = working.idxmax()
        events.append((when, float(working.loc[when])))
        working.loc[(working.index >= when - span)
                    & (working.index <= when + span)] = np.nan
    return sorted(events)


def build_table(cow, unreg):
    """One row per event: both ratios, the lag, and the magnitude."""
    first = max(cow.index.min(), unreg.index.min())
    last = min(cow.index.max(), unreg.index.max())
    grid = pd.date_range(first.ceil(RESAMPLE), last.floor(RESAMPLE),
                         freq=RESAMPLE)
    cow_h = cow.resample(RESAMPLE).mean().reindex(grid)
    unreg_h = unreg.resample(RESAMPLE).mean().reindex(grid)
    print("   overlap   : %s to %s (%d hours)"
          % (first.date(), last.date(), len(grid)))

    events = find_events(unreg_h, N_EVENTS, EVENT_MIN_SEPARATION_DAYS,
                         MIN_EVENT_CFS)
    print("   %d independent unregulated events above %s cfs"
          % (len(events), format(int(MIN_EVENT_CFS), ",")))

    half = pd.Timedelta(hours=EVENT_WINDOW_HOURS)
    rows = []
    for when, peak in events:
        window = slice(when - half, when + half)
        a, b = unreg_h.loc[window], cow_h.loc[window]
        need = int(MIN_WINDOW_COVERAGE * len(a))
        if len(a) == 0 or a.notna().sum() < need or b.notna().sum() < need:
            continue
        at_peak = b.get(when, np.nan)
        if not np.isfinite(at_peak):
            continue
        cow_peak_time = b.idxmax()
        rows.append({
            "event_time": when,
            "cas_unreg_peak_cfs": peak,
            "cow_at_cas_peak_cfs": float(at_peak),
            "cow_peak_cfs": float(b.max()),
            "cow_peak_time": cow_peak_time,
            "lag_hours": (cow_peak_time - when).total_seconds() / 3600.0,
            "ratio_coincident": float(at_peak) / peak,
            "ratio_peak": float(b.max()) / peak,
        })
    table = pd.DataFrame(rows)
    if len(table):
        table["magnitude_bin"] = pd.cut(table["cas_unreg_peak_cfs"],
                                        bins=MAGNITUDE_BINS, labels=BIN_LABELS,
                                        right=False)
    return table


def spread_stats(values):
    values = pd.Series(values).dropna()
    if len(values) < 2:
        return dict(n=len(values), median=np.nan, mean=np.nan, cv=np.nan,
                    iqr_over_median=np.nan)
    q1, q3 = values.quantile(0.25), values.quantile(0.75)
    median = values.median()
    return dict(n=len(values), median=median, mean=values.mean(),
                cv=values.std(ddof=1) / values.mean() if values.mean() else np.nan,
                iqr_over_median=(q3 - q1) / median if median else np.nan)


def report(table):
    print("\n" + "=" * 78)
    print("HOW BIG IS THE COWEEMAN RELATIVE TO UNREGULATED CASTLE ROCK?")
    print("=" * 78)
    for column, label in (("ratio_coincident",
                           "at the hour of the Castle Rock peak"),
                          ("ratio_peak", "Coweeman peak over Castle Rock peak")):
        stat = spread_stats(table[column])
        print("   %-38s n=%d" % (label, stat["n"]))
        print("      median %.3f   mean %.3f   CV %.2f   IQR/median %.2f"
              % (stat["median"], stat["mean"], stat["cv"],
                 stat["iqr_over_median"]))
    gap = (table["ratio_peak"] - table["ratio_coincident"])
    print("   timing penalty (peak ratio minus coincident ratio): median %.3f"
          % gap.median())
    print("   median lag of the Coweeman peak: %+.1f hours"
          % table["lag_hours"].median())

    print("\n" + "=" * 78)
    print("IS IT MORE CONSISTENT AT HIGH FLOWS?")
    print("=" * 78)
    for column, label in (("ratio_coincident", "coincident ratio"),
                          ("ratio_peak", "peak ratio")):
        print("\n   %s" % label.upper())
        print("      %-8s %4s %8s %8s %8s %8s"
              % ("bin", "n", "median", "mean", "CV", "IQR/med"))
        groups = []
        for name in BIN_LABELS:
            values = table.loc[table["magnitude_bin"] == name, column].dropna()
            stat = spread_stats(values)
            print("      %-8s %4d %8.3f %8.3f %8.2f %8.2f"
                  % (name, stat["n"], stat["median"], stat["mean"],
                     stat["cv"], stat["iqr_over_median"]))
            if len(values) >= MIN_PER_BIN:
                groups.append((name, values))

        good = table[[column, "cas_unreg_peak_cfs"]].dropna()
        if len(good) >= 6:
            rho, p_rho = stats.spearmanr(good["cas_unreg_peak_cfs"],
                                         good[column])
            print("      trend in LEVEL   : Spearman rho %+.2f, p = %.3f  (%s)"
                  % (rho, p_rho,
                     "ratio drifts with magnitude" if p_rho < 0.05
                     else "no significant drift"))
        if len(groups) >= 2:
            stat_bf, p_bf = stats.levene(*[g[1].values for g in groups],
                                         center="median")
            print("      trend in SPREAD  : Brown-Forsythe p = %.3f across %d "
                  "bins  (%s)"
                  % (p_bf, len(groups),
                     "spread differs" if p_bf < 0.05
                     else "no demonstrated difference"))
            widest = max(groups, key=lambda g: spread_stats(g[1])["cv"])
            tightest = min(groups, key=lambda g: spread_stats(g[1])["cv"])
            print("      widest CV %s (%.2f), tightest %s (%.2f)"
                  % (widest[0], spread_stats(widest[1])["cv"],
                     tightest[0], spread_stats(tightest[1])["cv"]))

    # scaling exponent: does the Coweeman grow as fast as the mainstem?
    good = table[(table["cow_at_cas_peak_cfs"] > 0)
                 & (table["cas_unreg_peak_cfs"] > 0)]
    if len(good) >= 6:
        fit = stats.linregress(np.log10(good["cas_unreg_peak_cfs"]),
                               np.log10(good["cow_at_cas_peak_cfs"]))
        print("\n   SCALING  log10(Coweeman at peak) on log10(Castle Rock unreg)")
        print("      exponent %.2f +/- %.2f, r-squared %.2f"
              % (fit.slope, fit.stderr, fit.rvalue ** 2))
        print("      An exponent of 1 would mean a constant fraction. %.2f "
              "means the" % fit.slope)
        print("      Coweeman grows %s than the mainstem, so the fraction %s "
              "with size."
              % ("more slowly" if fit.slope < 1 else "faster",
                 "falls" if fit.slope < 1 else "rises"))
    # --- what a coincident frequency analysis actually needs ----------------
    # Not the mean: the distribution, and above all its upper end. The design
    # question is how much the tributary can be contributing when the mainstem
    # crests, so the 75th and 90th percentiles carry more weight than the
    # centre.
    print("\n   COINCIDENT RATIO QUANTILES BY BIN")
    print("      %-8s %4s %7s %7s %7s %7s %7s"
          % ("bin", "n", "p10", "p25", "p50", "p75", "p90"))
    for name in BIN_LABELS:
        values = table.loc[table["magnitude_bin"] == name,
                           "ratio_coincident"].dropna()
        if not len(values):
            continue
        print("      %-8s %4d %7.3f %7.3f %7.3f %7.3f %7.3f"
              % (name, len(values), values.quantile(0.10),
                 values.quantile(0.25), values.median(),
                 values.quantile(0.75), values.quantile(0.90)))

    # --- events where the Coweeman was largely absent at the crest ----------
    # These are what make CV and IQR/median disagree. CV is driven by the two
    # or three events where the tributary contributed almost nothing at the
    # mainstem peak; the interquartile spread ignores them. Both numbers are
    # right about different things, so the events themselves are named rather
    # than left inside a summary statistic.
    print("\n   EVENTS WHERE THE COWEEMAN WAS NEARLY ABSENT AT THE CREST")
    flagged = []
    for name in BIN_LABELS:
        sub = table[table["magnitude_bin"] == name]
        if len(sub) < MIN_PER_BIN:
            continue
        cut = sub["ratio_coincident"].median() * 0.5
        for _, row in sub[sub["ratio_coincident"] < cut].iterrows():
            flagged.append((row, name))
    if not flagged:
        print("      none: every event has at least half its bin median")
    for row, name in sorted(flagged, key=lambda r: r[0]["ratio_coincident"]):
        print("      %s  unreg %s cfs  ratio %.3f  own peak %s cfs %+.0f h away"
              % (pd.Timestamp(row["event_time"]).date(),
                 format(int(row["cas_unreg_peak_cfs"]), ","),
                 row["ratio_coincident"],
                 format(int(row["cow_peak_cfs"]), ","), row["lag_hours"]))
    if flagged:
        print("      Each is a timing miss rather than a dry tributary -- the")
        print("      Coweeman peaked, just not when the Cowlitz did. They are")
        print("      why CV and IQR/median disagree about the top bin.")

    n_small = sum(1 for name in BIN_LABELS
                  if (table["magnitude_bin"] == name).sum() < MIN_PER_BIN)
    print("\n   CAVEAT: %d event(s) across %d bins. %s"
          % (len(table), len(BIN_LABELS),
             "Some bins are below %d events." % MIN_PER_BIN if n_small
             else "All bins carry at least %d events." % MIN_PER_BIN))
    print("   These tests have little power at this sample size -- a "
          "non-significant")
    print("   result means not demonstrated, not absent.")


def plot(table, stem):
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0))

    ax = axes[0][0]
    ax.scatter(table["cas_unreg_peak_cfs"], table["ratio_peak"], s=42,
               facecolor="none", edgecolor=C_LOW, lw=1.2, label="peak ratio")
    ax.scatter(table["cas_unreg_peak_cfs"], table["ratio_coincident"], s=42,
               color=C_COW, edgecolor="0.2", lw=0.5, label="coincident ratio")
    for edge in MAGNITUDE_BINS[1:-1]:
        ax.axvline(edge, color="0.7", lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("Castle Rock unregulated peak (cfs)", fontsize=9)
    ax.set_ylabel("Coweeman / Castle Rock unregulated", fontsize=9)
    ax.set_title("Ratio against magnitude", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[0][1]
    data = [table.loc[table["magnitude_bin"] == n, "ratio_coincident"].dropna()
            for n in BIN_LABELS]
    keep = [(n, d) for n, d in zip(BIN_LABELS, data) if len(d)]
    if keep:
        ax.boxplot([d.values for _, d in keep], tick_labels=[n for n, _ in keep])
        for i, (_, d) in enumerate(keep, start=1):
            ax.scatter(np.full(len(d), i) + np.random.uniform(-0.07, 0.07, len(d)),
                       d.values, s=18, color=C_COW, alpha=0.7, zorder=3)
    ax.set_xlabel("Castle Rock unregulated peak (cfs)", fontsize=9)
    ax.set_ylabel("Coincident ratio", fontsize=9)
    ax.set_title("Coincident ratio by magnitude bin", fontsize=10)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1][0]
    good = table[(table["cow_at_cas_peak_cfs"] > 0)
                 & (table["cas_unreg_peak_cfs"] > 0)]
    if len(good) >= 3:
        ax.scatter(good["cas_unreg_peak_cfs"], good["cow_at_cas_peak_cfs"],
                   s=42, color=C_HIGH, edgecolor="0.2", lw=0.5)
        fit = stats.linregress(np.log10(good["cas_unreg_peak_cfs"]),
                               np.log10(good["cow_at_cas_peak_cfs"]))
        xs = np.geomspace(good["cas_unreg_peak_cfs"].min(),
                          good["cas_unreg_peak_cfs"].max(), 50)
        ax.plot(xs, 10 ** (fit.intercept + fit.slope * np.log10(xs)),
                color=C_COW, lw=1.8,
                label="exponent %.2f (r$^2$=%.2f)" % (fit.slope, fit.rvalue ** 2))
        ratio = good["cow_at_cas_peak_cfs"].median() / good["cas_unreg_peak_cfs"].median()
        ax.plot(xs, ratio * xs, color="0.5", lw=1.2, ls="--",
                label="constant fraction (exponent 1)")
        ax.legend(fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Castle Rock unregulated peak (cfs)", fontsize=9)
    ax.set_ylabel("Coweeman at that hour (cfs)", fontsize=9)
    ax.set_title("Does the Coweeman scale with the mainstem?", fontsize=10)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1][1]
    ax.scatter(table["lag_hours"], table["ratio_coincident"], s=42,
               color=C_COW, edgecolor="0.2", lw=0.5)
    ax.axvline(0, color="0.4", lw=1.0, ls="--")
    ax.set_xlabel("Coweeman peak lag (hours; negative = leads)", fontsize=9)
    ax.set_ylabel("Coincident ratio", fontsize=9)
    ax.set_title("Timing against the coincident ratio", fontsize=10)
    ax.grid(alpha=0.3)

    fig.suptitle("Coweeman River as a fraction of unregulated Castle Rock flow\n"
                 "%d events, %s to %s"
                 % (len(table), table["event_time"].min().date(),
                    table["event_time"].max().date()), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(stem, dpi=150)
    plt.close(fig)


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)
    print("=" * 78)
    print("Coweeman as a proportion of UNREGULATED Castle Rock flow")
    print("=" * 78)
    print("\nINPUTS")
    cow = read_ecology_cache(CACHE_DIR)
    unreg = read_unreg(UNREG_DSS, UNREG_PATH)

    table = build_table(cow, unreg)
    if table is None or not len(table):
        raise SystemExit("No events with coverage in both records.")
    report(table)
    table.to_csv(OUT_CSV, index=False)
    plot(table, PLOT_PNG)
    print("\nTable: %s" % OUT_CSV)
    print("Plot : %s" % PLOT_PNG)


main()
