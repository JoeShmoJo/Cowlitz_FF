#BelowConfluence_FlowFrequency.py
# -*- coding: utf-8 -*-
"""
FINAL regulated peak flow-frequency curve for the Cowlitz BELOW the Coweeman
confluence, adding the three ungaged-at-the-confluence locals: Coweeman,
Arkansas Creek and Ostrander Creek.

THE METHOD, IN ONE LINE

    combined(AEP) = CAS_reg(AEP)
                    + CAS_unreg(AEP) x LOCAL_RATIO x lag_factor(AEP)

    LOCAL_RATIO   = UNIT_RUNOFF_UPLIFT x (sum local DA / CAS DA)
    lag_factor    = 0.80 -> 0.50, smooth in z, centred on 1% AEP

WHY THE LOCAL SCALES OFF THE *UNREGULATED* CURVE
    The Coweeman does not know Riffe and Mossyrock exist. Its flow responds to
    the storm, and the storm is what the UNREGULATED Castle Rock curve
    describes. Scaling the local off CAS_reg instead would shrink the
    tributary in proportion to a reservoir upstream of it -- about 18% low at
    the 1,000-yr.

    Pairing CAS_unreg(AEP) with CAS_reg(AEP) is not a coincidence assumption:
    it is the same river in the same event, and reg was DERIVED from unreg by
    routing in #Unreg_Reg_Curve.py. Nothing here assumes the Coweeman is
    simultaneously at its own 1,000-yr event -- that assumption, which the
    earlier same-AEP methods carried, is gone.

WHERE THE TWO FACTORS COME FROM

    LOCAL_RATIO -- drainage-area ratio, lightly corrected.
        Measured over 76 storm events, WY2007-WY2019, the Coweeman's own peak
        runs at 1.11x its plain drainage-area share of the unregulated Castle
        Rock peak in the >60k bin (0.0590 measured vs 0.0532 by area). That
        1.11 is UNIT_RUNOFF_UPLIFT. Set it to 1.00 for a pure area ratio; the
        two differ by ~10%, which is inside the scatter.

        Arkansas and Ostrander have no gage, so they inherit the Coweeman's
        unit rate. They are adjacent basins in the same low coastal hills, so
        that transfer is defensible. It would NOT be defensible from a
        dissimilar basin: East Fork Lewis, in the Cascade foothills, yields
        about twice the Coweeman's rate per unit area.

        The measured ratio is a LOWER BOUND. Three events were dropped
        because the Coweeman was above its rating at the crest and Ecology
        reported nothing -- those are events where the tributary was
        LARGEST. See Modules/ecology_io.py.

    lag_factor -- the tributary is not at its own crest when the mainstem is.
        The Coweeman peaks a median 13-20 h BEFORE the regulated Castle Rock
        crest and is already receding when the mainstem arrives. Measured
        ratio of Coweeman flow at the regulated crest to its own event peak:

            20-40k   n=51   median 0.809
            40-60k   n=18   median 0.762
            >60k     n= 7   median 0.420   (mean 0.511, range 0.36-0.81)

        The decline with magnitude is real (rho -0.32, p=0.004) but the tail
        value is poorly constrained: n=7, and the mechanism is a longer lead
        time (20 h vs 13 h) which is NOT itself significantly related to event
        size (rho -0.07, p=0.53). So 0.50 is a defensible tail estimate, not a
        measured constant. LAG_RARE is the single number to move; the
        sensitivity table at the foot of the run shows what it costs.

INPUT
    ../output/regulated_frequency_inferred.csv   (from #Unreg_Reg_Curve.py)
        supplies CAS_unreg and CAS_reg on one AEP grid, with 95% bounds.

OUTPUT
    ../output/below_confluence_frequency.csv
    ../output/diagnostics/below_confluence_frequency.png
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
OUT_CSV = r"../output/below_confluence_frequency.csv"
PLOT_PNG = r"../output/diagnostics/below_confluence_frequency.png"

CAS_DA = 2238.0                  # Cowlitz above Castle Rock, sq mi
TRIB_DA = {                      # sq mi at the confluence, not at the gage
    "Coweeman": 127.0,           # gage is 119; CDID3 used the same 1.07 step
    "Arkansas": 44.7,
    "Ostrander": 25.8,
}

UNIT_RUNOFF_UPLIFT = 1.11        # 1.00 = plain drainage-area ratio.
                                 # 1.11 = measured Coweeman rate, >60k bin.

LAG_COMMON = 0.80                # lag factor at common AEPs
LAG_RARE = 0.50                  # lag factor at rare AEPs -- least certain
                                 # number in this script, n=7. See docstring.
LAG_CENTER_AEP = 0.01            # logistic midpoint
LAG_WIDTH_Z = 0.40               # smaller = sharper transition, in z units

TARGET_AEP = 0.001               # called out on the plot and in the print
LAG_SENSITIVITY = (0.36, 0.50, 0.66, 0.80)   # tail range actually observed

C_REG = "#1a4f8a"
C_COMB = "#b7410e"
C_LOCAL = "#4c8c4a"

# ----------------------------------------------------------------------------


def lag_factor(aep):
    """Smooth LAG_COMMON -> LAG_RARE as the event gets rarer.

    Logistic in z rather than a step, so the combined curve has no kink at the
    transition. At LAG_CENTER_AEP it sits exactly halfway between the two.
    """
    z = stats.norm.ppf(1 - np.asarray(aep, dtype=float))
    z0 = stats.norm.ppf(1 - LAG_CENTER_AEP)
    return LAG_RARE + (LAG_COMMON - LAG_RARE) / (1 + np.exp((z - z0) / LAG_WIDTH_Z))


def local_ratio():
    """Local peak as a fraction of the unregulated Castle Rock peak."""
    return {name: UNIT_RUNOFF_UPLIFT * da / CAS_DA for name, da in TRIB_DA.items()}


def build(reg):
    ratios = local_ratio()
    out = pd.DataFrame({"AEP": reg["AEP"]})
    out["lag_factor"] = lag_factor(reg["AEP"].values)
    out["cowlitz_unreg_cfs"] = reg["unreg_computed_cfs"]
    out["cowlitz_reg_cfs"] = reg["reg_inferred_cfs"]

    local_total = np.zeros(len(reg))
    for name, ratio in ratios.items():
        peak = reg["unreg_computed_cfs"].values * ratio
        at_crest = peak * out["lag_factor"].values
        out["%s_peak_cfs" % name.lower()] = peak
        out["%s_at_crest_cfs" % name.lower()] = at_crest
        local_total += at_crest
    out["local_total_cfs"] = local_total
    out["combined_cfs"] = out["cowlitz_reg_cfs"] + out["local_total_cfs"]

    # Band: the regulated bound plus the local computed from the CORRESPONDING
    # unregulated bound, so a single hydrologic state drives both terms rather
    # than mixing a low mainstem with a high tributary.
    share = sum(ratios.values())
    for side, regcol, unregcol in [("lower", "reg_lower_95pct_cfs",
                                    "unreg_lower_95pct_cfs"),
                                   ("upper", "reg_upper_95pct_cfs",
                                    "unreg_upper_95pct_cfs")]:
        out["combined_%s_cfs" % side] = (
            reg[regcol].values
            + reg[unregcol].values * share * out["lag_factor"].values)

    out["local_pct_of_combined"] = 100 * out["local_total_cfs"] / out["combined_cfs"]
    out["increase_over_cowlitz_pct"] = (
        100 * out["local_total_cfs"] / out["cowlitz_reg_cfs"])
    return out, ratios


def report(out, ratios):
    share = sum(ratios.values())
    print("LOCAL RATIO  (local peak / unregulated Castle Rock peak)")
    for name, r in sorted(ratios.items(), key=lambda kv: -kv[1]):
        print("   %-10s %6.1f sq mi   %.5f   (%.4f by area x %.2f uplift)"
              % (name, TRIB_DA[name], r, TRIB_DA[name] / CAS_DA,
                 UNIT_RUNOFF_UPLIFT))
    print("   %-10s %6.1f sq mi   %.5f" % ("TOTAL", sum(TRIB_DA.values()), share))
    print("\n   lag factor: %.2f at common AEPs -> %.2f at rare, centred on "
          "%.2f%% AEP" % (LAG_COMMON, LAG_RARE, LAG_CENTER_AEP * 100))

    print("\n%8s %7s %11s %11s %10s %11s %7s"
          % ("AEP", "lag", "CAS unreg", "CAS reg", "local", "combined", "+%"))
    for _, r in out.iterrows():
        print("%8.4f %7.3f %11s %11s %10s %11s %7.1f"
              % (r["AEP"], r["lag_factor"],
                 format(int(r["cowlitz_unreg_cfs"]), ","),
                 format(int(r["cowlitz_reg_cfs"]), ","),
                 format(int(r["local_total_cfs"]), ","),
                 format(int(r["combined_cfs"]), ","),
                 r["increase_over_cowlitz_pct"]))

    row = out.iloc[(out["AEP"] - TARGET_AEP).abs().idxmin()]
    print("\nAt AEP=%.4f (1-in-%d):" % (row["AEP"], round(1 / row["AEP"])))
    print("   Cowlitz regulated      %11s cfs" % format(int(row["cowlitz_reg_cfs"]), ","))
    print("   local at the crest     %11s cfs  (%.1f%% of combined)"
          % (format(int(row["local_total_cfs"]), ","), row["local_pct_of_combined"]))
    print("   COMBINED               %11s cfs" % format(int(row["combined_cfs"]), ","))
    print("   95%% band               %11s to %s cfs"
          % (format(int(row["combined_lower_cfs"]), ","),
             format(int(row["combined_upper_cfs"]), ",")))

    print("\nSENSITIVITY to LAG_RARE at AEP=%.4f  (band width is %s cfs)"
          % (TARGET_AEP,
             format(int(row["combined_upper_cfs"] - row["combined_lower_cfs"]), ",")))
    base = row["combined_cfs"]
    keep = LAG_RARE
    for trial in LAG_SENSITIVITY:
        globals()["LAG_RARE"] = trial
        f = float(lag_factor(row["AEP"]))
        combined = row["cowlitz_reg_cfs"] + row["cowlitz_unreg_cfs"] * sum(
            ratios.values()) * f
        print("   LAG_RARE %.2f -> lag %.3f -> combined %11s cfs  (%+.1f%% vs base)"
              % (trial, f, format(int(combined), ","), 100 * (combined - base) / base))
    globals()["LAG_RARE"] = keep


def plot(out):
    z = stats.norm.ppf(1 - out["AEP"].values)
    fig, (ax, axl) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.4, 1]))

    ax.set_yscale("log")   # BEFORE any annotation -- a pre-log get_ylim()
                           # grabs a linear autoscale limit and crushes the
                           # data into a sliver once the scale changes.
    ax.fill_between(z, out["combined_lower_cfs"], out["combined_upper_cfs"],
                    color=C_COMB, alpha=0.13, label="Combined 95% band")
    ax.plot(z, out["cowlitz_unreg_cfs"], color="#9bb8d4", lw=1.4, ls=":",
            label="Cowlitz unregulated (drives the local)")
    ax.plot(z, out["cowlitz_reg_cfs"], color=C_REG, lw=2, ls="--",
            label="Cowlitz regulated at Castle Rock")
    ax.plot(z, out["combined_cfs"], color=C_COMB, lw=2.6,
            label="Below Coweeman confluence (regulated)")
    ax.axvline(stats.norm.ppf(1 - TARGET_AEP), color="gray", lw=1, ls=":")
    ax.text(stats.norm.ppf(1 - TARGET_AEP), 0.02, " 1,000-yr",
            transform=ax.get_xaxis_transform(), rotation=90, va="bottom",
            ha="right", fontsize=8, color="gray")
    ax.set_ylabel("Peak flow (cfs)")
    ax.set_title("Regulated peak flow frequency below the Coweeman confluence\n"
                 "Cowlitz regulated + Coweeman + Arkansas + Ostrander,\n"
                 "drainage-area scaled off the unregulated curve, then lagged",
                 fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    ticks = [0.99, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(stats.norm.ppf(1 - np.array(ticks)))
    top.set_xticklabels(["%g%%" % (a * 100) for a in ticks], rotation=45, fontsize=8)
    top.set_xlabel("AEP")

    bottom = np.zeros(len(out))
    for name, color in [("Coweeman", C_LOCAL), ("Arkansas", "#7fb069"),
                        ("Ostrander", "#b5d99c")]:
        v = out["%s_at_crest_cfs" % name.lower()].values
        axl.fill_between(z, bottom, bottom + v, color=color, alpha=0.85,
                         label="%s (%.1f sq mi)" % (name, TRIB_DA[name]))
        bottom = bottom + v
    axl.set_ylabel("Local at the crest (cfs)")
    axl.set_xlabel("Standard normal variate  z = $\\Phi^{-1}$(1 − AEP)")
    axl.grid(True, alpha=0.3)
    axl.legend(loc="upper left", fontsize=8, framealpha=0.95)

    axf = axl.twinx()
    axf.plot(z, out["lag_factor"], color="0.35", lw=1.6, ls=":")
    axf.set_ylim(0, 1)
    axf.set_ylabel("lag factor (dotted)", color="0.35", fontsize=9)

    fig.tight_layout()
    os.makedirs(os.path.dirname(PLOT_PNG), exist_ok=True)
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", PLOT_PNG)


def main():
    reg = pd.read_csv(REG_CSV).sort_values("AEP", ascending=False).reset_index(drop=True)
    out, ratios = build(reg)
    out.to_csv(OUT_CSV, index=False)
    report(out, ratios)
    plot(out)
    print("Wrote", OUT_CSV)


if __name__ == "__main__":
    main()
