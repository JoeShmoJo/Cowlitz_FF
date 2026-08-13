#Synthetic_Diagnostics.py
# -*- coding: utf-8 -*-
"""
Assess the synthetic ensemble design after the ResSim results are extracted.

The synthetics are a 4 x 4 x 3 factorial -- source event (shape), target
magnitude, and starting pool basis. This script answers the questions that
decide whether the design was the right one, and whether any axis can be
dropped or needs more members.

    1. VERIFICATION   Did the scaled unregulated peaks land on their targets?
                      A miss means the scale factor was applied to a different
                      quantity than the one being compared.

    2. WHICH AXIS MATTERS
                      At a fixed magnitude, how much does the regulated peak
                      move with SHAPE versus with STARTING POOL? Reported as
                      the spread in each direction. If pool spread is near zero
                      the pool axis is wasted; if shape spread is large the
                      curve cannot be a single line and needs a band.

    3. CURVE PLACEMENT
                      Synthetics plotted on the unreg-reg curve alongside the
                      period-of-record points, so it is visible whether they
                      extend the observed trend or depart from it.

    4. CONTROL LOSS   The local inflow is a hard floor on the regulated peak --
                      the project cannot release less than what enters below
                      the dam. Members approaching that floor are the ones
                      where the reservoir has run out of influence, and where
                      the curve should flatten.

    5. SCALING STRAIN Scale factors far from 1 assume hydrograph shape is
                      preserved with magnitude. Members are flagged by factor
                      so that assumption is visible rather than buried.

Run AFTER #Extract_Ensemble_To_Timeseries.py with SET_NAME = "ResSim_Synth".
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
from datetime import datetime
from pydsstools.heclib.dss import HecDss

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
SYNTH_DSS = r"../output/ResSim_Synth.dss"
SYNTH_MAPPING = r"../output/ensemble_synthetic_mapping.csv"
PATH_REG = "//CASTLEROCK_NWS/FLOW/*/1Hour/ResSim_Synth/"
PATH_LOCAL = "//CASTLEROCK_NWS/FLOW-LOCAL/*/1Hour/ResSim_Synth/"
PATH_MOS_IN = "//MOSSYROCK-POOL/FLOW-IN/*/1Hour/ResSim_Synth/"
PATH_POOL = "//MOSSYROCK-POOL/ELEV/*/1Hour/ResSim_Synth/"   # optional
# The unregulated flow AT CASTLE ROCK, routed by ResSim. This is the correct
# partner for the regulated peak: both are the same quantity at the same point,
# so their ratio is the reservoir effect and nothing else. Reported everywhere,
# and it is what #Unreg_Reg_Curve.py pairs with reg_peak.
PATH_UNREG = "//CASTLEROCK_NWS/FLOW-UNREG/*/1Hour/ResSim_Synth/"

# ALIGNMENT CHECK BASIS -- deliberately NOT the routed record.
# The member was scaled so that Mossyrock inflow + Castle Rock local, summed
# hour by hour and UNROUTED, hits the target exactly. Routing then attenuates
# that by 1-3%, varying by event and by magnitude. Checking the ROUTED peak
# against the target would be measuring routing, and it would need a tolerance
# loose enough to let a genuinely mis-read member slip through. The unrouted sum
# reproduces the target to 0.0%, which makes this an identity test: any
# departure at all means a wrong member was read.
# So: Flow-UNREG is the number reported; the unrouted sum is the number checked.
ALIGNMENT_TOLERANCE_FRAC = 0.005

# Period-of-record context, for the curve placement plot
POR_DATASET = r"../output/diagnostics/critical_duration_adjusted_dataset.csv"
POR_FIT_CSV = r"../output/critical_duration_adjusted_fits.csv"

OUT_CSV = r"../output/synthetic_results.csv"
DIAG_DIR = r"../output/diagnostics"
PLOT_STEM = r"../output/diagnostics/synthetic"

MAX_POOL_ELEV = 778.5
TARGET_TOLERANCE = 0.02      # fractional miss on the unregulated peak worth flagging
SCALE_STRAIN_WARN = 1.5      # scale factors beyond this are called out

# ----------------------------------------------------------------------------


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


def try_read(dss_file, pathname):
    try:
        return read_dss_series(dss_file, pathname)
    except Exception:
        return None


def block_stats(series, start, end):
    """Peak and its time within one member's synthetic-year block."""
    if series is None:
        return np.nan, pd.NaT
    window = series.loc[start:end].dropna()
    if len(window) == 0:
        return np.nan, pd.NaT
    return float(window.max()), pd.Timestamp(window.idxmax())


def spread_table(results):
    """At each magnitude, how far does the regulated peak move along each axis?"""
    rows = []
    for target, group in results.groupby("target", sort=False):
        good = group.dropna(subset=["reg_peak"])
        if len(good) < 2:
            continue
        # spread attributable to shape: range of per-event means
        by_event = good.groupby("event")["reg_peak"].mean()
        # spread attributable to pool: mean within-event range across pools
        by_pool = good.groupby("event")["reg_peak"].agg(lambda v: v.max() - v.min())
        rows.append({
            "target": target,
            "target_unreg_peak": good["target_unreg_peak_cfs"].iloc[0],
            "n": len(good),
            "reg_peak_mean": good["reg_peak"].mean(),
            "shape_spread_cfs": float(by_event.max() - by_event.min()),
            "pool_spread_cfs": float(by_pool.mean()),
            "total_spread_cfs": float(good["reg_peak"].max() - good["reg_peak"].min()),
        })
    table = pd.DataFrame(rows)
    if len(table):
        table["shape_pct"] = 100 * table["shape_spread_cfs"] / table["reg_peak_mean"]
        table["pool_pct"] = 100 * table["pool_spread_cfs"] / table["reg_peak_mean"]
        table["shape_over_pool"] = table["shape_spread_cfs"] / table["pool_spread_cfs"].replace(0, np.nan)
    return table


def plot_design(results, spreads, stem):
    """Regulated peak by magnitude, split out by event and by pool basis."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    targets = list(dict.fromkeys(results["target"]))
    xpos = {t: i for i, t in enumerate(targets)}

    ax = axes[0]
    for event, group in results.groupby("event"):
        g = group.dropna(subset=["reg_peak"]).groupby("target")["reg_peak"].mean()
        ax.plot([xpos[t] for t in g.index], g.values, marker="o", ms=6, lw=1.4,
                label=event)
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets)
    ax.set_ylabel("Regulated peak at Castle Rock (cfs)")
    ax.set_title("By source event (shape)\nmean across pool bases")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for basis, group in results.groupby("pool_basis"):
        g = group.dropna(subset=["reg_peak"]).groupby("target")["reg_peak"].mean()
        ax.plot([xpos[t] for t in g.index], g.values, marker="s", ms=6, lw=1.4,
                label=basis)
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets)
    ax.set_ylabel("Regulated peak at Castle Rock (cfs)")
    ax.set_title("By starting pool basis\nmean across source events")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.suptitle("Which design axis moves the regulated peak?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("%s_design.png" % stem, dpi=150)
    plt.close(fig)


def plot_curve(results, stem):
    """Synthetics on the unreg-reg curve with the period-of-record points."""
    fig, ax = plt.subplots(figsize=(10.5, 8))
    if os.path.isfile(POR_DATASET):
        por = pd.read_csv(POR_DATASET)
        col = "unreg_Peak_1-hr"
        if col in por and "reg_peak" in por:
            ax.loglog(por[col], por["reg_peak"], ls="none", marker="o", ms=5,
                      color="0.55", mec="0.3", mew=0.4, label="Period of record")
    markers = {"rulecurve": "^", "duration50": "s", "observed": "o"}
    colors = {}
    palette = ["#c0392b", "#2c7fb8", "#16a085", "#e67e22"]
    for i, event in enumerate(dict.fromkeys(results["event"])):
        colors[event] = palette[i % len(palette)]
    for (event, basis), group in results.groupby(["event", "pool_basis"]):
        good = group.dropna(subset=["reg_peak"])
        if len(good) == 0:
            continue
        ax.loglog(good["scaled_unreg_peak_cfs"], good["reg_peak"], ls="none",
                  marker=markers.get(basis, "d"), ms=8, color=colors[event],
                  mec="k", mew=0.5, alpha=0.85)
    good = results.dropna(subset=["reg_peak", "local_peak"])
    if len(good):
        order = good.sort_values("scaled_unreg_peak_cfs")
        ax.loglog(order["scaled_unreg_peak_cfs"], order["local_peak"],
                  color="k", lw=1.2, ls=":", label="Local inflow floor")
    handles = [Line2D([], [], color=colors[e], marker="o", ls="none", label=e)
               for e in colors]
    handles += [Line2D([], [], color="0.3", marker=m, ls="none", label=b)
                for b, m in markers.items()]
    handles += [Line2D([], [], color="0.55", marker="o", ls="none",
                       label="Period of record"),
                Line2D([], [], color="k", lw=1.2, ls=":", label="Local inflow floor")]
    ax.legend(handles=handles, fontsize=8, loc="upper left", ncol=2)
    ax.set_xlabel("Unregulated peak at Castle Rock (cfs)")
    ax.set_ylabel("Regulated peak at Castle Rock (cfs)")
    ax.set_title("Unregulated-regulated curve: synthetics against the record\n"
                 "colour = source event, symbol = starting pool basis")
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig("%s_curve.png" % stem, dpi=150)
    plt.close(fig)


def check_member_alignment(results):
    """Abort if a member's routed unregulated peak is not the peak it was built to.

    THE CHECK THAT MATTERS. Every synthetic member was scaled to hit an exact
    unregulated peak, and ResSim routes that same inflow, so the unregulated
    peak read back out has to match the target to within rounding. When it does
    not, the numbers are being read from the WRONG MEMBER, and the regulated
    peaks -- the whole point of the exercise -- belong to a different flood.

    That is not hypothetical. The synthetic results committed on 13 Aug 2026
    missed their targets by 3% to 29% and were used anyway; reading the members
    straight out of CAS_Synthetics.dss showed the regulated peaks were wrong by
    up to 130,000 cfs and were non-monotonic in every source event.

    The usual cause is stale members. ResSim does not clear old ensemble members
    when a smaller ensemble is run into the same simulation, so member numbers
    from a previous, larger run survive and the extract happily reads them.
    Rebuild the ResSim simulation, or delete the old members, and re-run.
    """
    if "target_miss_frac" not in results.columns:
        return
    # A member whose unregulated peak fell back to the target cannot be checked:
    # the comparison would be the target against itself. Treat that as a
    # failure, not a pass -- a check that cannot fail is not a check.
    if "unreg_peak_source" in results.columns:
        blind = results[results["unreg_peak_source"] != "sim"]
        if len(blind):
            print("=" * 78)
            print("MEMBER ALIGNMENT CANNOT BE CHECKED -- %d of %d members"
                  % (len(blind), len(results)))
            print("=" * 78)
            print("No simulated unregulated peak for these members, so there is")
            print("nothing independent to compare the target against.")
            print("\nThis needs Flow-Local at Castle Rock and Flow-IN at Mossyrock,")
            print("summed hour by hour -- NOT a Flow-UNREG record. Both are ordinary")
            print("ResSim outputs, so there is no need to make ResSim write all")
            print("records. Check that these two paths exist in the reassembled DSS:")
            print("   %s" % PATH_LOCAL)
            print("   %s" % PATH_MOS_IN)
            raise SystemExit("alignment check has no simulated peak to test against")

    miss = results["target_miss_frac"].abs()
    bad = results[miss > ALIGNMENT_TOLERANCE_FRAC]
    worst = float(miss.max()) if len(miss.dropna()) else np.nan
    if bad.empty:
        print("Member alignment  : OK -- every member's UNROUTED unregulated peak")
        print("                    matches its build target (worst miss %.2f%%, "
              "tolerance %.2f%%)"
              % (100 * worst, 100 * ALIGNMENT_TOLERANCE_FRAC))
        return
    print("=" * 78)
    print("MEMBER ALIGNMENT FAILED -- %d of %d members" % (len(bad), len(results)))
    print("=" * 78)
    print("Each member was BUILT to an exact unregulated peak. ResSim routes that")
    print("same inflow, so the unregulated peak read back must match. It does not,")
    print("which means these results were read from the wrong ensemble members and")
    print("the regulated peaks belong to different floods.\n")
    for _, row in bad.iterrows():
        print("   member %-3d %-9s %-7s target %10s  read back %10s   %+.1f%%"
              % (row["member"], row["event"], row["target"],
                 format(int(row["target_unreg_peak_cfs"]), ","),
                 format(int(row["unreg_peak_from_sim"]), ","),
                 100 * row["target_miss_frac"]))
    print("\nMost likely cause: STALE MEMBERS. ResSim does not clear old ensemble")
    print("members when a smaller ensemble is run into the same simulation, so")
    print("member numbers from a previous run survive alongside the new ones and")
    print("the extract reads whatever is at that number.")
    print("\nTo fix: confirm the member count in the ResSim simulation matches the")
    print("mapping CSV (%d members), clear the old members, re-run ResSim, then")
    print("re-run #Extract_Ensemble_To_Timeseries.py and this script.")
    print("\n%s was still written, for inspection. Do NOT use it." % OUT_CSV)
    raise SystemExit("member alignment check failed -- results not usable")


def main():
    if not os.path.isdir(DIAG_DIR):
        os.makedirs(DIAG_DIR)
    if not os.path.isfile(SYNTH_DSS):
        print("ERROR: %s not found." % SYNTH_DSS)
        print("Run #Extract_Ensemble_To_Timeseries.py with SET_NAME='ResSim_Synth' first.")
        return

    mapping = pd.read_csv(SYNTH_MAPPING, parse_dates=["real_start", "real_end",
                                                      "source_start"])
    reg = read_dss_series(SYNTH_DSS, PATH_REG)
    local = try_read(SYNTH_DSS, PATH_LOCAL)
    mos_in = try_read(SYNTH_DSS, PATH_MOS_IN)
    pool = try_read(SYNTH_DSS, PATH_POOL)
    unreg_routed = try_read(SYNTH_DSS, PATH_UNREG)

    rows = []
    for _, row in mapping.iterrows():
        start, end = row["real_start"], row["real_end"]
        reg_peak, reg_time = block_stats(reg, start, end)
        loc_peak, _ = block_stats(local, start, end)
        mos_peak, _ = block_stats(mos_in, start, end)
        pool_max, pool_time = block_stats(pool, start, end)
        # Sum the two series HOUR BY HOUR before taking the max. Adding the
        # independent maxima overstates the peak whenever the reservoir inflow
        # and the local do not crest in the same hour.
        unreg_peak, unreg_time = np.nan, pd.NaT
        if local is not None and mos_in is not None:
            a = local.loc[start:end]
            b = mos_in.loc[start:end]
            joint = (a.reindex(a.index.union(b.index)).fillna(0.0)
                     + b.reindex(a.index.union(b.index)).fillna(0.0)).dropna()
            if len(joint):
                unreg_peak, unreg_time = float(joint.max()), pd.Timestamp(joint.idxmax())
        # If the unregulated peak cannot be derived from the SIMULATION, the
        # alignment check has nothing independent to test against. Falling back
        # to the target here would make the check compare the target with
        # itself and pass every time, no matter which members were read -- so
        # record the fallback instead of hiding it.
        routed_peak, routed_time = block_stats(unreg_routed, start, end)
        unreg_source = "sim"
        if not np.isfinite(unreg_peak):
            unreg_peak = row["scaled_unreg_peak_cfs"]
            unreg_source = "target-fallback"
        entry = dict(row)
        entry.update({
            "reg_peak": reg_peak, "reg_peak_time": reg_time,
            "local_peak": loc_peak, "mos_in_peak": mos_peak,
            "unreg_peak_unrouted_cfs": unreg_peak,
            "unreg_peak_source": unreg_source,
            "unreg_peak_time": unreg_time,
            "unreg_peak_routed_cfs": routed_peak,
            "unreg_peak_routed_time": routed_time,
            "unreg_peak_from_sim": routed_peak if np.isfinite(routed_peak)
            else unreg_peak,
            "routing_loss_pct": (100.0 * (1.0 - routed_peak / unreg_peak)
                                 if np.isfinite(routed_peak) and unreg_peak
                                 else np.nan),
            "atten_ratio_routed": (reg_peak / routed_peak
                                   if np.isfinite(reg_peak) and routed_peak
                                   else np.nan),
            "sum_of_maxima_cfs": loc_peak + mos_peak
            if np.isfinite(loc_peak) and np.isfinite(mos_peak) else np.nan,
            # unrouted on purpose -- see ALIGNMENT_TOLERANCE_FRAC
            "target_miss_frac": (unreg_peak / row["target_unreg_peak_cfs"] - 1.0)
            if np.isfinite(unreg_peak) else np.nan,
            "atten_ratio": reg_peak / unreg_peak
            if np.isfinite(reg_peak) and unreg_peak else np.nan,
            "reg_over_local": reg_peak / loc_peak
            if np.isfinite(reg_peak) and np.isfinite(loc_peak) and loc_peak else np.nan,
            "max_pool_ft": pool_max,
            "pool_hit_max": bool(np.isfinite(pool_max) and pool_max >= MAX_POOL_ELEV),
        })
        rows.append(entry)
    results = pd.DataFrame(rows)
    results.to_csv(OUT_CSV, index=False, float_format="%.2f")

    check_member_alignment(results)
    if "routing_loss_pct" in results.columns and results["routing_loss_pct"].notna().any():
        loss = results["routing_loss_pct"]
        print("Routing to Castle Rock: unregulated peak attenuates %.2f%% to %.2f%%"
              % (loss.min(), loss.max()))
        print("                    (reported pairs use the routed peak; the")
        print("                     alignment check uses the unrouted sum)")

    got = results.dropna(subset=["reg_peak"])
    print("=" * 78)
    print("Members in mapping : %d      with results: %d" % (len(results), len(got)))
    if len(got) == 0:
        print("No results read -- check PATH_REG and the extraction step.")
        return

    print("\n1. TARGET VERIFICATION")
    miss = results.dropna(subset=["target_miss_frac"])
    bad = miss[miss["target_miss_frac"].abs() > TARGET_TOLERANCE]
    print("   scaled unregulated peak vs target: median %+.2f%%, worst %+.2f%%"
          % (100 * miss["target_miss_frac"].median(),
             100 * miss["target_miss_frac"].abs().max()))
    if len(bad):
        print("   %d members miss by more than %.0f%%:" % (len(bad), 100 * TARGET_TOLERANCE))
        for _, r in bad.head(6).iterrows():
            print("      %-8s %-7s %-11s %+.1f%%" % (r["event"], r["target"],
                                                     r["pool_basis"],
                                                     100 * r["target_miss_frac"]))

    print("\n2. WHICH AXIS MOVES THE ANSWER")
    spreads = spread_table(got)
    if len(spreads):
        print(spreads[["target", "n", "reg_peak_mean", "shape_spread_cfs",
                       "pool_spread_cfs", "shape_pct", "pool_pct",
                       "shape_over_pool"]].round(1).to_string(index=False))
        med = spreads["shape_over_pool"].median()
        if np.isfinite(med):
            if med > 2:
                print("   Shape dominates (%.1fx pool). The curve needs a band, not a"
                      % med)
                print("   line, and more shape variety would pay off more than more pools.")
            elif med < 0.5:
                print("   Starting pool dominates (%.1fx shape). The pool axis is doing"
                      % (1 / med))
                print("   the work -- consider a Monte Carlo over seasonal pool instead.")
            else:
                print("   Shape and pool contribute comparably (%.1fx). Keep both axes."
                      % med)
        spreads.to_csv(os.path.join(DIAG_DIR, "synthetic_axis_spread.csv"),
                       index=False, float_format="%.2f")

    print("\n3. CONTROL LOSS (regulated peak against the local inflow floor)")
    near = got.dropna(subset=["reg_over_local"])
    if len(near):
        print("   reg/local: median %.2f, min %.2f  (1.00 = the project is releasing"
              % (near["reg_over_local"].median(), near["reg_over_local"].min()))
        print("   nothing at the peak and the local alone sets the flow)")
        floored = near[near["reg_over_local"] < 1.05]
        print("   members within 5%% of the floor: %d of %d" % (len(floored), len(near)))
        for _, r in floored.head(6).iterrows():
            print("      %-8s %-7s %-11s reg %.0f vs local %.0f"
                  % (r["event"], r["target"], r["pool_basis"], r["reg_peak"],
                     r["local_peak"]))
    if got["pool_hit_max"].any():
        n = int(got["pool_hit_max"].sum())
        print("   %d members fill to max pool (%.1f ft) -- flood storage exhausted"
              % (n, MAX_POOL_ELEV))

    print("\n4. SCALING STRAIN")
    strained = results[results["scale_factor"] > SCALE_STRAIN_WARN]
    print("   scale factors: median %.2f, max %.2f" % (results["scale_factor"].median(),
                                                       results["scale_factor"].max()))
    print("   members scaled beyond %.1fx: %d (%s)"
          % (SCALE_STRAIN_WARN, len(strained),
             ", ".join(sorted(set(strained["event"]))) if len(strained) else "none"))

    plot_design(got, spreads, PLOT_STEM)
    plot_curve(got, PLOT_STEM)
    print("-" * 78)
    print("Results CSV : %s" % OUT_CSV)
    print("Spread CSV  : %s/synthetic_axis_spread.csv" % DIAG_DIR)
    print("Plots       : %s_design.png, %s_curve.png" % (PLOT_STEM, PLOT_STEM))


main()
