#Coincident_PerfectCorrelation.py
# -*- coding: utf-8 -*-
"""
Castle Rock combined peak flow: Cowlitz (regulated) + Coweeman, assuming
perfect correlation.

WHAT THIS IS
    The simplest of the coincident-frequency options, and the one with the
    most direct precedent: CDID3's own Phase 1 report (USACE 2013a) used
    this exact method for the Coweeman/Cowlitz/Columbia combination --
    "perfect correlation was assumed between the Columbia, Cowlitz, and
    Coweeman Rivers. This was a conservative assumption to identify whether
    the levee would qualify under the most extreme conditions." (See
    CAS_Reg_Unreg/docs/CDID3_Coincident_Frequency_Notes.md.)

    Concretely: pair the SAME AEP from each independent curve and add the
    two flows. If the Cowlitz's 1% event is 118,000 cfs and the Coweeman's
    1% event is 11,000 cfs, the combined "1% ACE" flow is reported as
    129,000 cfs. There is no probability combination step -- the AEP label
    on the combined curve is inherited from the AEP the two inputs share,
    which is why this is a conservative bound rather than a real joint
    probability statement: actually getting BOTH rivers at their 1% event
    at the same time is rarer than 1% unless they are perfectly correlated,
    which they are not (see #Coincident_CorrConditioned.py for the
    correlation-conditioned alternative, and coweeman_proportion.csv for
    the actual measured correlation).

    This is NOT the same thing as the Oakridge-study envelope method,
    which pairs DIFFERENT AEPs across variables for an inundation
    footprint. Here every pair shares one AEP, and that AEP is asserted
    (not derived) to be the combined event's AEP.

WHY BUILD IT AT ALL, GIVEN IT'S CONSERVATIVE BY CONSTRUCTION
    It is real USACE precedent for this exact river system, it takes one
    formula to build once both marginal curves exist, and it gives a
    conservative bookend to compare the correlation-conditioned curve
    against: if the two methods come out close together, the correlation
    story isn't doing much work; if they're far apart, that gap IS the
    story of how much crediting realistic (imperfect) correlation actually
    saves.

UNCERTAINTY
    Perfect correlation (r=1) between two random variables A and B means
    SD(A+B) = SD(A) + SD(B) -- NOT the root-sum-of-squares combination
    used elsewhere in this repo for INDEPENDENT uncertainty terms (see
    #Unreg_Reg_Curve.py). Summing the lower bounds and summing the upper
    bounds arithmetically is the internally consistent way to propagate
    uncertainty under this method's own r=1 assumption, so that's what
    this script does -- not because the two curves' actual uncertainty is
    perfectly correlated (it may not be) but because that is what "assume
    perfect correlation" means once you also want an uncertainty band, and
    this script should not smuggle in a different assumption for the band
    than it uses for the central estimate.

INPUTS
    CAS_Reg_Unreg/output/regulated_frequency_inferred.csv
        reg_inferred_cfs, reg_lower_95pct_cfs, reg_upper_95pct_cfs
        (from #Unreg_Reg_Curve.py)
    CAS_Reg_Unreg/output/diagnostics/coweeman_frequency_table.csv
        Value, LowerConf, UpperConf
        (from #Coweeman_FlowFrequency.py -- see that script's docstring for
        why this curve is an interim approximation, not an HEC-SSP result)
    Both are already on the same 16-point AEP grid (0.99 down to 0.0001),
    so this script does no interpolation -- it is a row-for-row sum.
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
REG_CSV = r"../output/regulated_frequency_inferred.csv"
COW_CSV = r"../output/diagnostics/coweeman_frequency_table.csv"

OUT_DIR = r"../output/diagnostics"
OUT_CSV = os.path.join(OUT_DIR, "coincident_perfect_correlation.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coincident_perfect_correlation.png")

TARGET_AEP = 0.001     # 1,000-year -- the study's stated target, called out
                        # on the plot; both curves already cover down to
                        # 0.0001 (10,000-year) so there is margin either way.

C_REG = "#1a4f8a"
C_COMBINED = "#b7410e"
C_COW = "#4c8c4a"

# ----------------------------------------------------------------------------


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    reg = pd.read_csv(REG_CSV)
    cow = pd.read_csv(COW_CSV)

    if not np.allclose(sorted(reg["AEP"]), sorted(cow["AEP"])):
        raise ValueError(
            "AEP grids differ between %s and %s -- this script sums row for "
            "row and assumes they already match. Interpolate one onto the "
            "other before combining rather than silently misaligning rows."
            % (REG_CSV, COW_CSV))

    reg = reg.sort_values("AEP", ascending=False).reset_index(drop=True)
    cow = cow.sort_values("AEP", ascending=False).reset_index(drop=True)
    assert (reg["AEP"].values == cow["AEP"].values).all()

    out = pd.DataFrame({"AEP": reg["AEP"]})
    out["cowlitz_reg_cfs"] = reg["reg_inferred_cfs"]
    out["cowlitz_reg_lower_cfs"] = reg["reg_lower_95pct_cfs"]
    out["cowlitz_reg_upper_cfs"] = reg["reg_upper_95pct_cfs"]
    out["coweeman_cfs"] = cow["Value"]
    out["coweeman_lower_cfs"] = cow["LowerConf"]
    out["coweeman_upper_cfs"] = cow["UpperConf"]
    out["combined_cfs"] = out["cowlitz_reg_cfs"] + out["coweeman_cfs"]
    # Perfect-correlation (r=1) propagation: SD(A+B) = SD(A) + SD(B), i.e.
    # sum the bounds directly rather than root-sum-of-squares. See the
    # UNCERTAINTY section of the module docstring.
    out["combined_lower_cfs"] = out["cowlitz_reg_lower_cfs"] + out["coweeman_lower_cfs"]
    out["combined_upper_cfs"] = out["cowlitz_reg_upper_cfs"] + out["coweeman_upper_cfs"]
    out["coweeman_pct_of_combined"] = 100 * out["coweeman_cfs"] / out["combined_cfs"]

    out.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV)
    print(out[["AEP", "cowlitz_reg_cfs", "coweeman_cfs", "combined_cfs",
               "coweeman_pct_of_combined"]].to_string(index=False))

    target_row = out.iloc[(out["AEP"] - TARGET_AEP).abs().idxmin()]
    print("\nAt AEP=%.4f (%s): Cowlitz reg=%.0f cfs + Coweeman=%.0f cfs = "
          "combined %.0f cfs (Coweeman is %.1f%% of the total)"
          % (target_row["AEP"], "1,000-yr" if TARGET_AEP == 0.001 else "target",
             target_row["cowlitz_reg_cfs"], target_row["coweeman_cfs"],
             target_row["combined_cfs"], target_row["coweeman_pct_of_combined"]))

    # -- plot --
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(z, out["cowlitz_reg_cfs"], color=C_REG, lw=2, ls="--",
            label="Cowlitz regulated alone (Castle Rock)")
    ax.plot(z, out["coweeman_cfs"], color=C_COW, lw=1.5, ls=":",
            label="Coweeman alone (interim curve)")
    ax.plot(z, out["combined_cfs"], color=C_COMBINED, lw=2.5,
            label="Combined, same-AEP sum (perfect-correlation assumption)")
    ax.fill_between(z, out["combined_lower_cfs"], out["combined_upper_cfs"],
                     color=C_COMBINED, alpha=0.15,
                     label="Combined 5-95% band (bounds summed, r=1)")
    ax.set_yscale("log")
    ax.axvline(stats.norm.ppf(1 - TARGET_AEP), color="gray", lw=1, ls=":")
    # Anchored in axes-fraction y so this doesn't depend on ylim being
    # finalized yet -- see #Coincident_CorrConditioned.py for why placing it
    # via a pre-log-scale get_ylim() crushes the real data into a sliver.
    ax.text(stats.norm.ppf(1 - TARGET_AEP), 0.01, " 1,000-yr", transform=ax.get_xaxis_transform(),
            rotation=90, va="bottom", ha="right", fontsize=8, color="gray")
    ax.set_xlabel("Standard normal variate  (z = Φ⁻¹(1 − AEP))")
    ax.set_ylabel("Flow at Castle Rock confluence (cfs)")
    ax.set_title("Coincident Castle Rock peak — perfect-correlation (same-AEP) method\n"
                  "CDID3 Phase 1 (2013) precedent — conservative bound, not a joint probability")
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
