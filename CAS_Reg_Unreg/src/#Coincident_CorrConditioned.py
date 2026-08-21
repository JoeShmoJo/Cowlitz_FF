#Coincident_CorrConditioned.py
# -*- coding: utf-8 -*-
"""
Castle Rock combined peak flow: correlation-conditioned combination.

WHAT THIS IS
    The classical middle path between "assume perfect correlation"
    (#Coincident_PerfectCorrelation.py, the CDID3 Phase 1 precedent) and a
    full CDID3-Phase-2a-style Monte Carlo (correlation matrix + timing
    shape sets + hydraulic routing -- more machinery than a flow-only
    combined peak needs, see CAS_Reg_Unreg/docs/
    CDID3_Coincident_Frequency_Notes.md).

    Both marginal curves (Cowlitz unregulated -> regulated, Coweeman) are
    built independently, same as the perfect-correlation method. What
    differs is the pairing: instead of summing the SAME AEP from each
    curve, this estimates the correlation between the two rivers' PEAK
    magnitudes (in z-space, from the 79-event concurrent record already
    built by #Coweeman_Proportion.py) and uses it to find the
    CONDITIONALLY EXPECTED value of one river given the other's AEP, then
    sums. A Cowlitz 1% event isn't paired with the Coweeman's own 1% event
    (that overstates how often both are that extreme at once, unless
    r=1) -- it's paired with whatever Coweeman flow the correlation says
    is typical when the Cowlitz is at that magnitude.

    Checked in BOTH directions (Cowlitz-controls, Coweeman-controls) and
    enveloped, per the classical convention -- neither stream is assumed
    to always be the one driving the combined event.

WHY PEAK-TO-PEAK CORRELATION, NOT THE COINCIDENT-HOUR VALUE
    coweeman_proportion.csv carries two different pairings, deliberately:
    ratio_peak (each river's own peak, wherever in the event window it
    falls) and ratio_coincident (the Coweeman's flow at the EXACT hour the
    Cowlitz peaks, which already has the observed ~10-20 hour lag baked
    in). This script uses ratio_peak's underlying values (cow_peak_cfs vs
    cas_unreg_peak_cfs) -- the classical method's job is to capture
    magnitude correlation with a single, stable parameter that is safe to
    hold fixed out to a 1,000-year target; folding in the timing lag as
    well would mean extrapolating a second, thinner relationship (see the
    "regression on raw events" method this repo's chat history rejected
    for exactly that reason). This method is consequently ANOTHER
    perfect-timing assumption, same as the perfect-correlation method --
    it improves on that method's magnitude assumption (correlated instead
    of equal-AEP) but not its timing assumption. Worth stating plainly,
    not hiding.

METHOD
    1. For each of the 79 paired events, convert cas_unreg_peak_cfs to a
       z-statistic via the Cowlitz UNREGULATED curve, and cow_peak_cfs to a
       z-statistic via the Coweeman curve (both by log-flow-vs-z
       interpolation on each curve's own fitted points).
    2. r = Pearson correlation of the two z-series.
    3. Direction A (Cowlitz controls): for each target AEP, z_cas =
       Phi^-1(1-AEP); conditional z_cow = r * z_cas (bivariate-normal
       conditional mean for standardized marginals); invert through the
       Coweeman curve; combined = Cowlitz REGULATED flow at that AEP (via
       #Unreg_Reg_Curve.py's transform) + conditional Coweeman flow.
    4. Direction B (Coweeman controls): mirror image -- condition the
       Cowlitz's UNREGULATED flow on the Coweeman's AEP, push it through
       the unreg->reg transform, add the Coweeman's own flow at that AEP.
    5. Combined = max(Direction A, Direction B) at each AEP -- the
       classical envelope-both-directions step.

VALID RANGE -- FLOOD TAIL ONLY, SAME LIMIT CDID3 STATES FOR ITS OWN CURVES
    Conditioning on a correlation regresses the OTHER variable toward ITS
    OWN median whenever |r|<1 -- expected behavior, but it means that at a
    very common AEP (e.g. 0.99, a small, frequent event), the "conditionally
    expected" partner value sits closer to its typical flow than a naive
    same-percentile scaling would suggest. The visible symptom: at AEP=0.99
    this method's envelope comes out HIGHER than the perfect-correlation
    method's same-AEP sum -- the opposite of the flood-tail relationship,
    where this method is correctly LOWER (about 2% lower at both 1% and
    0.1% AEP in the current run) than assuming perfect correlation. That
    is not a bug in this script; it is the same thing CDID3 says about its
    own frequency curves -- "results should not be used to infer the water
    surface elevations for smaller events (e.g. 50% ACE)" -- extended to
    the combination step. Trust this method's output for flood-tail AEPs
    (roughly 10% and rarer); the crossover near the common end is a visible
    marker of where that stops being true, left in the plot rather than
    trimmed out.

CORRELATION SAMPLE SIZE -- SAME CAVEAT AS EVERYWHERE ELSE IN THIS ANALYSIS
    r is estimated from n=79 events, only 9 of which exceed 60,000 cfs at
    Castle Rock. It is held FIXED across the whole AEP grid, including out
    to 1,000-year and beyond -- that is a much smaller extrapolation than
    stretching an empirical regression curve's shape that far (the reason
    this method was preferred over that one), but it is still an
    extrapolation of a single parameter past where it was measured, and
    should be reported as such, not as settled.

INPUTS
    CAS_Reg_Unreg/output/diagnostics/coweeman_proportion.csv   (event pairs)
    CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv          (Cowlitz unreg)
    CAS_Reg_Unreg/output/regulated_frequency_inferred.csv      (Cowlitz reg
                                                                 + the unreg->
                                                                 reg pairing)
    CAS_Reg_Unreg/output/diagnostics/coweeman_frequency_table.csv (Coweeman)
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
PROPORTION_CSV = r"../output/diagnostics/coweeman_proportion.csv"
UNREG_CSV = r"../../CAS_Unreg_FF/output/CAS_Unreg_frequency_table.csv"
REG_CSV = r"../output/regulated_frequency_inferred.csv"
COW_CSV = r"../output/diagnostics/coweeman_frequency_table.csv"
PERFECT_CORR_CSV = r"../output/diagnostics/coincident_perfect_correlation.csv"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "coincident_corr_conditioned.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coincident_corr_conditioned.png")

TARGET_AEP = 0.001     # 1,000-year

C_ENV = "#b7410e"
C_DIRA = "#4c8c4a"
C_DIRB = "#8a5aa8"
C_PERFECT = "gray"

# ----------------------------------------------------------------------------


def z_from_aep(aep):
    return stats.norm.ppf(1 - np.asarray(aep, dtype=float))


def build_flow_z_curve(aep, flow):
    """Return (log10_flow, z) sorted ascending by z -- ready for np.interp
    in either direction, since both are monotonic in a real frequency curve."""
    z = z_from_aep(aep)
    order = np.argsort(z)
    return np.log10(np.asarray(flow)[order]), z[order]


def z_of_flow(flow_value, log_flow_grid, z_grid, extrap_count):
    lf = np.log10(flow_value)
    if lf < log_flow_grid[0] or lf > log_flow_grid[-1]:
        extrap_count[0] += 1
    return np.interp(lf, log_flow_grid, z_grid)


def flow_of_z(z_value, log_flow_grid, z_grid):
    return 10 ** np.interp(z_value, z_grid, log_flow_grid)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    events = pd.read_csv(PROPORTION_CSV)
    unreg = pd.read_csv(UNREG_CSV)
    unreg = unreg[unreg["Duration"] == "Peak"].copy()
    reg = pd.read_csv(REG_CSV)
    cow = pd.read_csv(COW_CSV)

    unreg_lf, unreg_z = build_flow_z_curve(unreg["AEP"], unreg["Value"])
    cow_lf, cow_z = build_flow_z_curve(cow["AEP"], cow["Value"])

    # -- 1-2: correlation of peak z-statistics over the concurrent record --
    n_extrap_cas, n_extrap_cow = [0], [0]
    z_cas = np.array([z_of_flow(v, unreg_lf, unreg_z, n_extrap_cas) for v in events["cas_unreg_peak_cfs"]])
    z_cow = np.array([z_of_flow(v, cow_lf, cow_z, n_extrap_cow) for v in events["cow_peak_cfs"]])
    r, pval = stats.pearsonr(z_cas, z_cow)
    n = len(events)
    n_large = int((events["cas_unreg_peak_cfs"] > 60000).sum())
    print("Peak-to-peak correlation (z-space): r=%.3f (p=%.4f), n=%d events "
          "(%d exceed 60,000 cfs at Castle Rock)" % (r, pval, n, n_large))
    if n_extrap_cow[0]:
        print("  %d/%d Coweeman peaks fall below the curve's 99%%-AEP point "
              "(2,023 cfs) -- most of these 79 events are common, not flood-"
              "tail, events, so this is expected linear extrapolation at the "
              "small end, not a curve-range problem." % (n_extrap_cow[0], n))
    if n_extrap_cas[0]:
        print("  %d/%d Cowlitz peaks required extrapolation past the "
              "unregulated curve's fitted range." % (n_extrap_cas[0], n))

    # -- unreg -> reg transform, from the paired columns #Unreg_Reg_Curve.py
    #    already wrote (monotonic increasing in both, log-log interpolation) --
    reg_sorted = reg.sort_values("unreg_computed_cfs")
    unreg_to_reg_lf = np.log10(reg_sorted["unreg_computed_cfs"].values)
    unreg_to_reg_reg = reg_sorted["reg_inferred_cfs"].values

    def reg_from_unreg(unreg_cfs):
        lf = np.log10(unreg_cfs)
        return np.interp(lf, unreg_to_reg_lf, np.log10(unreg_to_reg_reg))

    def reg_from_unreg_cfs(unreg_cfs):
        return 10 ** reg_from_unreg(unreg_cfs)

    # -- 3-5: both directions, on the shared AEP grid --
    aep_grid = sorted(reg["AEP"].values, reverse=True)
    rows = []
    for aep in aep_grid:
        z_target = z_from_aep(aep)

        # Direction A: Cowlitz controls.
        z_cow_cond = r * z_target
        coweeman_cond_cfs = flow_of_z(z_cow_cond, cow_lf, cow_z)
        cowlitz_reg_at_aep = float(reg.loc[np.isclose(reg["AEP"], aep), "reg_inferred_cfs"].iloc[0])
        combined_a = cowlitz_reg_at_aep + coweeman_cond_cfs

        # Direction B: Coweeman controls.
        z_cas_cond = r * z_target
        cas_unreg_cond_cfs = flow_of_z(z_cas_cond, unreg_lf, unreg_z)
        cas_reg_cond_cfs = reg_from_unreg_cfs(cas_unreg_cond_cfs)
        coweeman_at_aep = float(cow.loc[np.isclose(cow["AEP"], aep), "Value"].iloc[0])
        combined_b = cas_reg_cond_cfs + coweeman_at_aep

        combined_env = max(combined_a, combined_b)
        controlling = "Cowlitz" if combined_a >= combined_b else "Coweeman"

        rows.append({
            "AEP": aep,
            "z_target": z_target,
            "correlation_r": r,
            "dirA_cowlitz_reg_cfs": cowlitz_reg_at_aep,
            "dirA_coweeman_conditional_cfs": coweeman_cond_cfs,
            "dirA_combined_cfs": combined_a,
            "dirB_cowlitz_reg_conditional_cfs": cas_reg_cond_cfs,
            "dirB_coweeman_cfs": coweeman_at_aep,
            "dirB_combined_cfs": combined_b,
            "combined_envelope_cfs": combined_env,
            "controlling_direction": controlling,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV)
    print(out[["AEP", "dirA_combined_cfs", "dirB_combined_cfs",
               "combined_envelope_cfs", "controlling_direction"]].to_string(index=False))

    target_row = out.iloc[(out["AEP"] - TARGET_AEP).abs().idxmin()]
    print("\nAt AEP=%.4f (1,000-yr): envelope combined = %.0f cfs (%s-controlled)"
          % (target_row["AEP"], target_row["combined_envelope_cfs"],
             target_row["controlling_direction"]))

    # -- plot, with the perfect-correlation curve overlaid for comparison --
    z = z_from_aep(out["AEP"].values)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(z, out["dirA_combined_cfs"], color=C_DIRA, lw=1.5, ls="--",
            label="Direction A: Cowlitz controls")
    ax.plot(z, out["dirB_combined_cfs"], color=C_DIRB, lw=1.5, ls="--",
            label="Direction B: Coweeman controls")
    ax.plot(z, out["combined_envelope_cfs"], color=C_ENV, lw=2.5,
            label="Envelope (max of A, B) — correlation-conditioned, r=%.2f" % r)

    if os.path.exists(PERFECT_CORR_CSV):
        pc = pd.read_csv(PERFECT_CORR_CSV)
        ax.plot(z_from_aep(pc["AEP"].values), pc["combined_cfs"], color=C_PERFECT,
                lw=1.5, ls=":", label="Perfect-correlation method (for comparison)")

    ax.set_yscale("log")
    ax.axvline(z_from_aep(TARGET_AEP), color="gray", lw=1, ls=":")
    # Anchored in axes-fraction y (not data y) so this doesn't depend on
    # ylim being finalized yet -- placing it via a pre-log-scale get_ylim()
    # grabs a linear-scale (often negative) value that then forces the log
    # axis to stretch to contain it, crushing the real data into a sliver.
    ax.text(z_from_aep(TARGET_AEP), 0.01, " 1,000-yr", transform=ax.get_xaxis_transform(),
            rotation=90, va="bottom", ha="right", fontsize=8, color="gray")
    ax.set_xlabel("Standard normal variate  (z = Φ⁻¹(1 − AEP))")
    ax.set_ylabel("Flow at Castle Rock confluence (cfs)")
    ax.set_title("Coincident Castle Rock peak — correlation-conditioned method\n"
                  "r=%.2f from n=%d concurrent events (%d exceed 60,000 cfs)"
                  % (r, n, n_large))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    aep_ticks = [0.99, 0.5, 0.1, 0.02, 0.01, 0.002, 0.001, 0.0002, 0.0001]
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(z_from_aep(np.array(aep_ticks)))
    ax2.set_xticklabels(["%.2f%%" % (a * 100) for a in aep_ticks], rotation=45, fontsize=8)
    ax2.set_xlabel("AEP")

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("Wrote", PLOT_PNG)


if __name__ == "__main__":
    main()
