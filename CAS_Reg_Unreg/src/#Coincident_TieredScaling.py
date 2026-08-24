#Coincident_TieredScaling.py
# -*- coding: utf-8 -*-
"""
Castle Rock regulated peak + Arkansas Creek + Ostrander Creek + Coweeman,
combined with a two-tier same-AEP scaling factor.

THE RULE
    combined(AEP) = Cowlitz_reg(AEP) + factor(AEP) * tributary(AEP)
    Applied to each of Arkansas Creek, Ostrander Creek, and the Coweeman
    independently, then summed with Cowlitz's regulated curve. No separate
    timing-lag correction on top -- both endpoints of factor() are
    COINCIDENT ratios (tributary flow at the moment Castle Rock peaks,
    over the tributary's own peak), so whatever lag exists is already
    inside the number; adding a second timing adjustment would
    double-count it.

    factor(AEP) is a SMOOTH logistic transition in z-space (z =
    Phi^-1(1-AEP), the same axis every curve in this analysis is plotted
    against), not a hard step -- changed from an earlier version of this
    script after the step was flagged as an unrealistic discontinuity to
    write into a design curve. It runs from FACTOR_COMMON at common AEPs
    to FACTOR_RARE at rare AEPs, centered on TRANSITION_CENTER_AEP with
    width TRANSITION_WIDTH_Z (in z-units):

        factor(z) = FACTOR_RARE + (FACTOR_COMMON - FACTOR_RARE)
                    / (1 + exp((z - z_center) / TRANSITION_WIDTH_Z))

    TRANSITION_WIDTH_Z=0.4 puts most of the transition within about one
    AEP_GRID step on either side of the center (roughly the 2% to 0.5%
    AEP range) -- narrow enough to still mean something like "common" vs
    "rare" behavior, wide enough not to be a cliff. It is a judgment
    call, not a fitted number; there are only 3 magnitude bins per basin
    behind this, nowhere near enough to fit a transition shape to.

WHERE THE TWO NUMBERS COME FROM, AND HOW MUCH TO TRUST THEM
    0.80: the 2009-era East Fork Lewis analog-basin figure
    (COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx, Section B.6) -- median 0.809 in
    the Coweeman's OWN 20,000-40,000 cfs bin (n=51) independently matches
    this closely, so it's reasonably well supported in the moderate range
    across two different basins.

    0.50: the Coweeman's own >60,000 cfs bin (n=8, median 0.494-0.520
    depending on regulated- vs unregulated-timing basis -- see
    Downstream_Confluence_Notes.md sections 4-5). NOT measured at or
    beyond 1% AEP -- unregulated Castle Rock at 1% AEP is 168,884 cfs and
    this bin's largest event is 155,018 cfs. Using 0.50 for AEP<=1% is a
    deliberate extrapolation past the last actual data point, not an
    interpolation within it, chosen because a flat, simple, honestly-
    extrapolated number was judged more defensible than a finer curve
    fit to 8 points. Carried over to Arkansas and Ostrander with ZERO
    data of their own to check it against -- two stacked assumptions
    (the AEP extrapolation, and the cross-basin transfer), not one.

    UPDATE -- the East Fork Lewis check ran (#EFLewis_Analog_Check.py,
    n=71, real USGS data, no longer blocked once fetched locally). It
    corroborates the DIRECTION of the tail effect but not its MAGNITUDE:
        Coweeman:   20-40k 0.809 -> >60k 0.520  (n=8,  36% drop)
        EF Lewis:   20-40k 0.669 -> >60k 0.624  (n=18,  7% drop)
    EF Lewis's tail bin has more than double Coweeman's sample and comes
    from a USGS gage without Coweeman's documented rating-curve capping
    problem -- there is a real argument that 0.50 is too aggressive a
    reduction and something closer to 0.62-0.65 is both better evidenced
    and the more conservative choice (LOWER factor means LESS added
    flow, i.e. 0.50 is the less conservative of the two, not the
    safer one -- easy to get backwards). FACTOR_RARE is left at 0.50
    below pending an explicit decision on this; change the one constant
    when that's settled.

CAP AT AEP=0.001 (1,000-YEAR)
    Arkansas/Ostrander's StreamStats table (CDID3_Coincident... no --
    Downstream_Confluence_Notes.md Table B-VI) only extends to AEP=0.001.
    Cowlitz and Coweeman both support the full grid down to 0.0001, but
    this script does not extrapolate Arkansas/Ostrander past their source
    table, so the combined output stops at 0.001 -- which is exactly the
    1,000-year target this whole sub-task was scoped to, so nothing is
    lost for the stated purpose.

UNCERTAINTY -- HONEST GAP, NOT FILLED IN
    Cowlitz and Coweeman both carry a 5-95% band, propagated through the
    same tiered factor and summed. Arkansas and Ostrander have NO
    uncertainty information in their source table (StreamStats gave point
    estimates only) -- their contribution is added to the band's center
    only. The resulting combined band therefore UNDERSTATES true
    uncertainty; it is missing two of five contributing terms entirely.
    Not fixed here -- flagged in the output instead of quietly ignored.
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
OUT_CSV = os.path.join(OUT_DIR, "coincident_tiered_scaling.csv")
PLOT_PNG = os.path.join(OUT_DIR, "coincident_tiered_scaling.png")
SENS_CSV = os.path.join(OUT_DIR, "coincident_method_sensitivity.csv")

TRANSITION_CENTER_AEP = 0.01   # 1% -- logistic transition midpoint
TRANSITION_WIDTH_Z = 0.4       # smaller = sharper transition, in z-units
FACTOR_COMMON = 0.80           # factor at common AEPs (AEP >> center)
FACTOR_RARE = 0.50             # factor at rare AEPs (AEP << center) -- see
                                # docstring's EF Lewis update before trusting
                                # this over ~0.62-0.65

# Downstream_Confluence_Notes.md, Table B-VI (USGS StreamStats via
# COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx). Point estimates only -- no
# uncertainty band available for either creek.
ARKANSAS_CFS = {
    0.999: 410, 0.99: 580, 0.95: 790, 0.90: 930, 0.80: 1130, 0.70: 1300,
    0.60: 1470, 0.50: 1620, 0.40: 1800, 0.30: 2040, 0.20: 2300, 0.10: 2740,
    0.05: 3200, 0.02: 3770, 0.01: 4220, 0.005: 4650, 0.002: 5310, 0.001: 5700,
}
OSTRANDER_CFS = {
    0.999: 230, 0.99: 325, 0.95: 440, 0.90: 520, 0.80: 620, 0.70: 710,
    0.60: 790, 0.50: 881, 0.40: 960, 0.30: 1100, 0.20: 1220, 0.10: 1470,
    0.05: 1700, 0.02: 2010, 0.01: 2240, 0.005: 2480, 0.002: 2790, 0.001: 3010,
}
# Cap the combined output at the StreamStats table's own range.
MAX_RARE_AEP_KEY = min(ARKANSAS_CFS)   # 0.001

C_CAS = "#1a4f8a"
C_COMBINED = "#b7410e"
C_TIER = "gray"

# ----------------------------------------------------------------------------


def scaling_factor(aep):
    """Smooth logistic transition from FACTOR_COMMON to FACTOR_RARE as AEP
    gets rarer, centered (in z-space) on TRANSITION_CENTER_AEP. Replaces an
    earlier hard step at the same center -- see docstring."""
    z = stats.norm.ppf(1 - aep)
    z_center = stats.norm.ppf(1 - TRANSITION_CENTER_AEP)
    return FACTOR_RARE + (FACTOR_COMMON - FACTOR_RARE) / (
        1 + np.exp((z - z_center) / TRANSITION_WIDTH_Z))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    reg = pd.read_csv(REG_CSV)
    cow = pd.read_csv(COW_CSV)

    # Intersection of all three sources' AEP grids -- StreamStats carries
    # points (0.999, 0.70, 0.60, 0.40, 0.30) the reg/Coweeman curves don't
    # tabulate, and the reg/Coweeman curves extend past 0.001 to points
    # StreamStats doesn't cover. Only the shared points can be combined
    # without interpolating something.
    reg_aeps = set(np.round(reg["AEP"].values, 6))
    cow_aeps = set(np.round(cow["AEP"].values, 6))
    trib_aeps = set(ARKANSAS_CFS)
    aep_grid = sorted(reg_aeps & cow_aeps & trib_aeps, reverse=True)
    dropped_trib = sorted(trib_aeps - set(aep_grid), reverse=True)
    dropped_curve = sorted((reg_aeps & cow_aeps) - set(aep_grid), reverse=True)
    if dropped_trib:
        print("StreamStats AEPs not on the reg/Coweeman grid, skipped:", dropped_trib)
    if dropped_curve:
        print("reg/Coweeman AEPs beyond StreamStats' range, skipped:", dropped_curve)

    def lookup(df, aep, col):
        row = df.loc[np.isclose(df["AEP"], aep)]
        if len(row) == 0:
            raise ValueError("AEP %s not found in %s -- interpolation not implemented, "
                              "add the exact point or extend the source table instead."
                              % (aep, col))
        return float(row[col].iloc[0])

    rows = []
    for aep in aep_grid:
        factor = scaling_factor(aep)

        cas_best = lookup(reg, aep, "reg_inferred_cfs")
        cas_lo = lookup(reg, aep, "reg_lower_95pct_cfs")
        cas_hi = lookup(reg, aep, "reg_upper_95pct_cfs")

        ark = ARKANSAS_CFS[aep]
        ost = OSTRANDER_CFS[aep]

        cow_best = lookup(cow, aep, "Value")
        cow_lo = lookup(cow, aep, "LowerConf")
        cow_hi = lookup(cow, aep, "UpperConf")

        below_ark = cas_best + factor * ark
        below_ost = below_ark + factor * ost
        below_cow = below_ost + factor * cow_best

        # Band: Cowlitz's own band plus Coweeman's own band (arithmetic sum,
        # same "assume high correlation" logic as #Coincident_
        # PerfectCorrelation.py), scaled by the same tier factor. Arkansas/
        # Ostrander contribute to the center only -- see module docstring.
        band_lo = cas_lo + factor * (ark + ost + cow_lo)
        band_hi = cas_hi + factor * (ark + ost + cow_hi)

        rows.append({
            "AEP": aep,
            "scaling_factor": factor,
            "cowlitz_reg_cfs": cas_best,
            "arkansas_cfs_x_factor": factor * ark,
            "ostrander_cfs_x_factor": factor * ost,
            "coweeman_cfs_x_factor": factor * cow_best,
            "below_arkansas_cfs": below_ark,
            "below_ostrander_cfs": below_ost,
            "below_coweeman_cfs": below_cow,
            "combined_lower_incomplete_band_cfs": band_lo,
            "combined_upper_incomplete_band_cfs": band_hi,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV)
    print(out[["AEP", "scaling_factor", "cowlitz_reg_cfs", "below_arkansas_cfs",
               "below_ostrander_cfs", "below_coweeman_cfs"]].to_string(index=False))

    target = out.iloc[(out["AEP"] - 0.001).abs().idxmin()]
    print("\nAt AEP=%.4f (1,000-yr): Cowlitz reg=%.0f cfs -> below Coweeman confluence=%.0f cfs"
          % (target["AEP"], target["cowlitz_reg_cfs"], target["below_coweeman_cfs"]))
    print("(band shown in the CSV is INCOMPLETE -- missing Arkansas/Ostrander's own "
          "uncertainty entirely; see module docstring)")

    # -- plot: two panels --
    #   top    absolute curves, with each tributary's contribution shaded so
    #          Arkansas and Ostrander are actually visible rather than being
    #          a sliver hidden under the combined line
    #   bottom the same contributions as a PERCENT of the Cowlitz-only curve,
    #          which is the only way the two small creeks read at all on a
    #          log axis spanning 20,000-200,000 cfs
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, (ax, axp) = plt.subplots(2, 1, figsize=(9.5, 9), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})

    base = out["cowlitz_reg_cfs"].values
    a1 = out["below_arkansas_cfs"].values
    a2 = out["below_ostrander_cfs"].values
    a3 = out["below_coweeman_cfs"].values

    ax.fill_between(z, base, a1, color="#4c8c4a", alpha=0.55, label="+ Arkansas Creek")
    ax.fill_between(z, a1, a2, color="#8a5aa8", alpha=0.55, label="+ Ostrander Creek")
    ax.fill_between(z, a2, a3, color="#e08a4a", alpha=0.65, label="+ Coweeman River")
    ax.plot(z, base, color=C_CAS, lw=2, ls="--", label="Cowlitz regulated alone (Castle Rock)")
    ax.plot(z, a3, color=C_COMBINED, lw=2.5, label="Combined below Coweeman confluence")
    ax.fill_between(z, out["combined_lower_incomplete_band_cfs"],
                     out["combined_upper_incomplete_band_cfs"],
                     color=C_COMBINED, alpha=0.10,
                     label="Band (incomplete -- Cowlitz+Coweeman only)")
    ax.set_yscale("log")
    ax.set_ylabel("Flow below Coweeman confluence (cfs)")
    ax.set_title("Downstream combined peak -- smooth same-AEP scaling\n"
                  "factor(AEP): %.2f (common) -> %.2f (rare), logistic, centered "
                  "on %.0f%% AEP" % (FACTOR_COMMON, FACTOR_RARE, TRANSITION_CENTER_AEP * 100))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

    axp.fill_between(z, 0, 100 * (a1 - base) / base, color="#4c8c4a", alpha=0.55)
    axp.fill_between(z, 100 * (a1 - base) / base, 100 * (a2 - base) / base,
                      color="#8a5aa8", alpha=0.55)
    axp.fill_between(z, 100 * (a2 - base) / base, 100 * (a3 - base) / base,
                      color="#e08a4a", alpha=0.65)
    axp.plot(z, 100 * (a3 - base) / base, color=C_COMBINED, lw=2)
    axp.set_ylabel("Added flow, % of Cowlitz alone")
    axp.set_xlabel("Standard normal variate  (z = Phi^-1(1 - AEP))")
    axp.grid(True, alpha=0.3)

    # scaling_factor itself, on its own right-hand axis of the bottom panel,
    # so the smooth transition shape is directly visible, not just implied
    # by the curves it produced.
    axf = axp.twinx()
    axf.plot(z, out["scaling_factor"], color=C_TIER, lw=1.5, ls=":")
    axf.set_ylabel("scaling_factor(AEP)", color=C_TIER, fontsize=9)
    axf.tick_params(axis="y", colors=C_TIER, labelsize=8)
    axf.set_ylim(0, 1)

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(stats.norm.ppf(1 - np.array(aep_grid)))
    ax2.set_xticklabels(["%.2f%%" % (a * 100) for a in aep_grid], rotation=45, fontsize=7)
    ax2.set_xlabel("AEP")

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("Wrote", PLOT_PNG)

    # -- method sensitivity: does the factor choice actually matter? --
    # Recover each tributary's UNSCALED contribution by dividing out the tier
    # factor, then rebuild the two alternative methods from the same inputs.
    trib = (out["arkansas_cfs_x_factor"] + out["ostrander_cfs_x_factor"]
            + out["coweeman_cfs_x_factor"]) / out["scaling_factor"]
    m100 = out["cowlitz_reg_cfs"] + trib
    m80 = out["cowlitz_reg_cfs"] + 0.80 * trib
    mstep = out["below_coweeman_cfs"]
    band = (out["combined_upper_incomplete_band_cfs"]
            - out["combined_lower_incomplete_band_cfs"])
    sens = pd.DataFrame({
        "AEP": out["AEP"],
        "cowlitz_alone_cfs": out["cowlitz_reg_cfs"],
        "combined_stepped_cfs": mstep,
        "combined_80pct_cfs": m80,
        "combined_100pct_cfs": m100,
        "trib_share_of_100pct_pct": 100 * trib / m100,
        "method_spread_cfs": m100 - mstep,
        "uncertainty_band_cfs": band,
        "band_over_method_spread": band / (m100 - mstep),
    })
    sens.to_csv(SENS_CSV, index=False)
    print("Wrote", SENS_CSV)
    print()
    print("METHOD SENSITIVITY -- the factor choice vs. the uncertainty it sits inside:")
    print(sens[["AEP", "combined_stepped_cfs", "combined_80pct_cfs",
                "combined_100pct_cfs", "trib_share_of_100pct_pct",
                "method_spread_cfs", "band_over_method_spread"]]
          .round(1).to_string(index=False))
    tail = sens[np.isclose(sens["AEP"], 0.001)].iloc[0]
    print()
    print("At 1,000-yr: tributaries are %.1f%% of the combined total, but the whole "
          "spread between the\n  three scaling methods is %.0f cfs -- %.0fx SMALLER than "
          "the uncertainty band (%.0f cfs)\n  it sits inside. The factor choice is not "
          "what decides this number."
          % (tail["trib_share_of_100pct_pct"], tail["method_spread_cfs"],
             tail["band_over_method_spread"], tail["uncertainty_band_cfs"]))


if __name__ == "__main__":
    main()
