#Coincident_TieredScaling.py
# -*- coding: utf-8 -*-
"""
Castle Rock regulated peak + Arkansas Creek + Ostrander Creek + Coweeman,
combined with a two-tier same-AEP scaling factor.

THE RULE (decided in chat, not derived fresh here)
    For AEP > 1% (more frequent than 1-in-100):
        combined(AEP) = Cowlitz_reg(AEP) + 0.80 * tributary(AEP)
    For AEP <= 1% (1-in-100 and rarer, out to 1,000-yr):
        combined(AEP) = Cowlitz_reg(AEP) + 0.50 * tributary(AEP)
    Applied to each of Arkansas Creek, Ostrander Creek, and the Coweeman
    independently, then summed with Cowlitz's regulated curve. No separate
    timing-lag correction on top -- both factors are COINCIDENT ratios
    (tributary flow at the moment Castle Rock peaks, over the tributary's
    own peak), so whatever lag exists is already inside the number; adding
    a second timing adjustment would double-count it.

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

    A same-basin, better-instrumented check (East Fork Lewis itself,
    against the long unregulated Castle Rock series back to 1928) was
    proposed to see whether the same high-AEP/low-AEP transition shows
    up independently -- blocked in this environment (USGS hosts are not
    reachable from this sandbox's network policy). Worth running locally,
    not done here.

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

TIER_BREAK_AEP = 0.01     # 1% -- at and rarer than this, use FACTOR_RARE
FACTOR_COMMON = 0.80      # AEP > 0.01
FACTOR_RARE = 0.50        # AEP <= 0.01

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


def tier_factor(aep):
    return FACTOR_RARE if aep <= TIER_BREAK_AEP else FACTOR_COMMON


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
        factor = tier_factor(aep)

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
            "tier_factor": factor,
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
    print(out[["AEP", "tier_factor", "cowlitz_reg_cfs", "below_arkansas_cfs",
               "below_ostrander_cfs", "below_coweeman_cfs"]].to_string(index=False))

    target = out.iloc[(out["AEP"] - 0.001).abs().idxmin()]
    print("\nAt AEP=%.4f (1,000-yr): Cowlitz reg=%.0f cfs -> below Coweeman confluence=%.0f cfs"
          % (target["AEP"], target["cowlitz_reg_cfs"], target["below_coweeman_cfs"]))
    print("(band shown in the CSV is INCOMPLETE -- missing Arkansas/Ostrander's own "
          "uncertainty entirely; see module docstring)")

    # -- plot --
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(z, out["cowlitz_reg_cfs"], color=C_CAS, lw=2, ls="--",
            label="Cowlitz regulated alone (Castle Rock)")
    ax.plot(z, out["below_coweeman_cfs"], color=C_COMBINED, lw=2.5,
            label="Combined below Coweeman confluence (tiered 80%/50%)")
    ax.fill_between(z, out["combined_lower_incomplete_band_cfs"],
                     out["combined_upper_incomplete_band_cfs"],
                     color=C_COMBINED, alpha=0.15,
                     label="Band (incomplete -- Cowlitz+Coweeman only)")
    ax.axvline(stats.norm.ppf(1 - TIER_BREAK_AEP), color=C_TIER, lw=1, ls=":")
    ax.text(stats.norm.ppf(1 - TIER_BREAK_AEP), 0.02, " 1% tier break\n (80% -> 50%)",
            transform=ax.get_xaxis_transform(), rotation=90, va="bottom", ha="right",
            fontsize=8, color=C_TIER)

    ax.set_yscale("log")
    ax.set_xlabel("Standard normal variate  (z = Φ⁻¹(1 − AEP))")
    ax.set_ylabel("Flow below Coweeman confluence (cfs)")
    ax.set_title("Downstream combined peak — tiered same-AEP scaling\n"
                  "80% common (AEP>1%) / 50% rare (AEP<=1%), Arkansas+Ostrander+Coweeman")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    aep_ticks = aep_grid
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(stats.norm.ppf(1 - np.array(aep_ticks)))
    ax2.set_xticklabels(["%.2f%%" % (a * 100) for a in aep_ticks], rotation=45, fontsize=7)
    ax2.set_xlabel("AEP")

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("Wrote", PLOT_PNG)


if __name__ == "__main__":
    main()
