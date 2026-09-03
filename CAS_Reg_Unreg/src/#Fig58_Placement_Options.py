"""
Three ways to place the adjusted regulated peaks on Figure 5-8, side by side.

Written for the DQC comment C18 after the within window ranking made the
unregulated curve look like a poor fit. Each panel shows the same curves and
bands, the full 95 year unregulated record ranked as HEC-SSP draws it, and the
41 usable adjusted regulated peaks placed one of three ways:

  A  ranked among themselves over the 51 assessed years, WY1974 to WY2024.
     The literal reading of the reviewer's request. The cloud sits above the
     curve because the window holds nine of the ten largest unregulated years.
  B  ranked among themselves, then given the plotting positions their own 41
     water years hold in the 95 year unregulated record, largest to largest.
     Monotone, on the curve's own probability axis, centred on the curve. The
     rank preserving assumption is the one the transform itself makes.
  C  each peak at the plotting position its own water year holds in the 95
     year record, paired by year. Not monotone. The version the reviewer
     objected to.

Diagnostic only. Writes output/diagnostics/fig58_placement_options.png and
prints the median of observed over curve for each panel. Does not touch
Figure 5-8, which #Unreg_Reg_Curve.py produces.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

os.chdir(os.path.dirname(os.path.abspath(__file__)))

UNREG_RECORD_CSV = r"../../CAS_Unreg_FF/output/wy_record_ssp.csv"
ADJUSTED_PEAKS_CSV = r"../output/adjusted_peaks.csv"
REG_FREQ_CSV = r"../output/regulated_frequency_inferred.csv"
OUT_PNG = r"../output/diagnostics/fig58_placement_options.png"

WINDOW = (1974, 2024)
AEP_TICKS = [0.99, 0.9, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
FLOW_LIMITS = (10000.0, 400000.0)
C_UNREG, C_REG = "#2c7fb8", "#c0392b"


def median_pp(n_rank, n_record):
    i = np.arange(1, n_rank + 1)
    return (i - 0.3) / (n_record + 0.4)


def curve_at(freq, col, aep):
    z = norm.ppf(1.0 - freq["AEP"].values)
    y = np.log10(freq[col].values)
    o = np.argsort(z)
    return 10.0 ** np.interp(norm.ppf(1.0 - np.asarray(aep, dtype=float)), z[o], y[o])


def main():
    unreg = pd.read_csv(UNREG_RECORD_CSV)
    adj = pd.read_csv(ADJUSTED_PEAKS_CSV)
    freq = pd.read_csv(REG_FREQ_CSV).sort_values("AEP")
    usable = adj[adj.screen_passed == True].copy()
    n_assessed = len(adj)

    full = unreg.sort_values("Peak", ascending=False).reset_index(drop=True)
    full["pp"] = median_pp(len(full), len(full))
    pp_of_wy = dict(zip(full.WY, full.pp))

    reg_sorted = np.sort(usable.adjusted_peak.values)[::-1]
    own_pp = usable.WY.map(pp_of_wy).values

    panels = {
        "A  ranked within WY%d to WY%d (n = %d of %d)" % (WINDOW + (len(usable), n_assessed)):
            (median_pp(len(usable), n_assessed), reg_sorted),
        "B  ranked, at the 95 year positions of the same %d water years" % len(usable):
            (np.sort(own_pp), reg_sorted),
        "C  each peak at its own water year's 95 year position (current, not monotone)":
            (own_pp, usable.adjusted_peak.values),
    }

    zx = lambda p: norm.ppf(1.0 - np.asarray(p, dtype=float))
    a = freq["AEP"].values
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.6), sharey=True)
    for ax, (title, (pp, val)) in zip(axes, panels.items()):
        ratio = val / curve_at(freq, "reg_inferred_cfs", pp)
        print("%s: median obs/curve %.3f, %d of %d above" % (title, np.median(ratio), int((ratio > 1).sum()), len(val)))
        ax.fill_between(zx(a), freq["unreg_lower_90pct_cfs"], freq["unreg_upper_90pct_cfs"], color=C_UNREG, alpha=0.12, lw=0)
        ax.fill_between(zx(a), freq["reg_lower_90pct_cfs"], freq["reg_upper_90pct_cfs"], color=C_REG, alpha=0.12, lw=0)
        ax.plot(zx(a), freq["unreg_expected_cfs"], color=C_UNREG, lw=2.2, label="Unregulated")
        ax.plot(zx(a), freq["reg_inferred_cfs"], color=C_REG, lw=2.2, label="Regulated")
        ax.scatter(zx(full.pp), full.Peak, s=16, facecolors="none", edgecolors=C_UNREG, lw=0.9,
                   label="Unregulated record (n=95)")
        ax.scatter(zx(pp), val, s=20, facecolors="none", edgecolors=C_REG, lw=1.1,
                   label="Adjusted regulated peaks (n=%d)" % len(val))
        ax.set_yscale("log")
        ax.set_ylim(*FLOW_LIMITS)
        ax.set_xlim(zx(0.99), zx(0.001))
        ax.set_xticks(zx(AEP_TICKS))
        ax.set_xticklabels(["%g" % (100 * p) for p in AEP_TICKS])
        ax.set_xlabel("Annual exceedance probability (%)")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_title("%s\nmedian obs/curve %.2f, %d of %d above" % (title, np.median(ratio), int((ratio > 1).sum()), len(val)), fontsize=9.5)
    axes[0].set_ylabel("Peak flow (cfs)")
    axes[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "{:,.0f}".format(v)))
    axes[0].legend(loc="upper left", fontsize=8.5, framealpha=0.95)
    fig.suptitle("Figure 5-8 candidates. Same curves and bands, same unregulated record, three placements of the adjusted regulated peaks", fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=140)
    print("wrote", OUT_PNG)


if __name__ == "__main__":
    main()
