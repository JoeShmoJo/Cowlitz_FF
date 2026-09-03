"""
Diagnostic for the DQC comment on Figure 5-8 (memo comment C18).

The reviewer asked for the adjusted regulated peaks to be ranked among
themselves and plotted at plotting positions, the standard convention. When
that was done the cloud sat above the regulated curve, roughly 33 of 41 points
over the line, and looked as though it could not have produced the curve.

This figure shows why, using the unregulated record as the control. Both
records are ranked INSIDE the same 51 water year window, WY1974 to WY2024, at
median plotting positions (i - 0.3) / (n + 0.4), and drawn against the adopted
curves. The unregulated points from that window sit above the unregulated
curve by the same margin the regulated points sit above the regulated curve,
because the window holds 9 of the 10 largest unregulated years in the 95 year
record. A ranked cloud from a flood rich window plots above ANY curve fitted to
the full record. The offset is the sample window, not the transform.

Reads only published outputs. Writes one PNG to output/diagnostics and prints
the numbers quoted in the caption. Does not touch Figure 5-8 itself, which is
still produced by #Unreg_Reg_Curve.py.
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
OUT_PNG = r"../output/diagnostics/fig58_window_ranked_check.png"

WINDOW = (1974, 2024)          # the assessed regulated record
AEP_TICKS = [0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
FLOW_LIMITS = (10000.0, 400000.0)
C_UNREG, C_REG, C_FULL = "#2c7fb8", "#c0392b", "#9e9e9e"


def median_pp(values, n_record):
    """Median plotting positions, largest first, ranked over n_record slots."""
    v = np.sort(np.asarray(values, dtype=float))[::-1]
    i = np.arange(1, len(v) + 1)
    return (i - 0.3) / (n_record + 0.4), v


def curve_at(freq, col, aep):
    """Log flow interpolated on the normal variate of AEP."""
    z = norm.ppf(1.0 - freq["AEP"].values)
    y = np.log10(freq[col].values)
    o = np.argsort(z)
    return 10.0 ** np.interp(norm.ppf(1.0 - np.asarray(aep, dtype=float)), z[o], y[o])


def main():
    unreg = pd.read_csv(UNREG_RECORD_CSV)
    adj = pd.read_csv(ADJUSTED_PEAKS_CSV)
    freq = pd.read_csv(REG_FREQ_CSV).sort_values("AEP")

    win = unreg[(unreg.WY >= WINDOW[0]) & (unreg.WY <= WINDOW[1])]
    n_win = len(win)
    usable = adj[adj.screen_passed == True]
    n_assessed = len(adj)

    # ranked inside the window, both records
    u_aep, u_val = median_pp(win.Peak.values, n_win)
    r_aep, r_val = median_pp(usable.adjusted_peak.values, n_assessed)
    # full record, for contrast
    f_aep, f_val = median_pp(unreg.Peak.values, len(unreg))

    u_ratio = u_val / curve_at(freq, "unreg_expected_cfs", u_aep)
    r_ratio = r_val / curve_at(freq, "reg_inferred_cfs", r_aep)
    f_ratio = f_val / curve_at(freq, "unreg_expected_cfs", f_aep)

    top10 = unreg.nlargest(10, "Peak").WY
    in_win = int(((top10 >= WINDOW[0]) & (top10 <= WINDOW[1])).sum())

    print("Ranked inside WY%d-%d, median plotting positions" % WINDOW)
    print("  unregulated, n=%d over %d: median obs/curve %.3f, %d of %d above the curve"
          % (n_win, n_win, np.median(u_ratio), int((u_ratio > 1).sum()), n_win))
    print("  adjusted regulated, n=%d over %d: median obs/curve %.3f, %d of %d above the curve"
          % (len(usable), n_assessed, np.median(r_ratio), int((r_ratio > 1).sum()), len(usable)))
    print("Full unregulated record, n=%d: median obs/curve %.3f, %d of %d above"
          % (len(unreg), np.median(f_ratio), int((f_ratio > 1).sum()), len(unreg)))
    print("Largest ten unregulated years inside the window: %d of 10" % in_win)
    print("Window median unregulated peak %.0f cfs against %.0f cfs for the full record"
          % (win.Peak.median(), unreg.Peak.median()))

    # ---- figure
    fig, ax = plt.subplots(figsize=(11, 8))
    zx = lambda p: norm.ppf(1.0 - np.asarray(p, dtype=float))
    a = freq["AEP"].values
    ax.fill_between(zx(a), freq["unreg_lower_90pct_cfs"], freq["unreg_upper_90pct_cfs"],
                    color=C_UNREG, alpha=0.12, lw=0, label="Unregulated, 90% (HEC-SSP)")
    ax.fill_between(zx(a), freq["reg_lower_90pct_cfs"], freq["reg_upper_90pct_cfs"],
                    color=C_REG, alpha=0.12, lw=0, label="Regulated, 90% (frequency + transform)")
    ax.plot(zx(a), freq["unreg_expected_cfs"], color=C_UNREG, lw=2.5, label="Unregulated curve, fitted to WY1927-2026 (n=95)")
    ax.plot(zx(a), freq["reg_inferred_cfs"], color=C_REG, lw=2.5, label="Regulated curve, transform of the unregulated curve")
    ax.scatter(zx(f_aep), f_val, s=14, marker="x", color=C_FULL, lw=0.8,
               label="Unregulated record ranked over all 95 years")
    ax.scatter(zx(u_aep), u_val, s=26, facecolors="none", edgecolors=C_UNREG, lw=1.2,
               label="Unregulated WY%d-%d ranked within the window (n=%d)" % (WINDOW + (n_win,)))
    ax.scatter(zx(r_aep), r_val, s=26, facecolors="none", edgecolors=C_REG, lw=1.2,
               label="Adjusted regulated WY%d-%d ranked within the window (n=%d of %d)"
               % (WINDOW + (len(usable), n_assessed)))

    ax.set_yscale("log")
    ax.set_ylim(*FLOW_LIMITS)
    ax.set_xlim(zx(0.99), zx(0.001))
    ax.set_xticks(zx(AEP_TICKS))
    ax.set_xticklabels(["%g" % (100 * p) for p in AEP_TICKS])
    ax.set_xlabel("Annual exceedance probability (%)")
    ax.set_ylabel("Peak flow (cfs)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "{:,.0f}".format(v)))
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title("Both records ranked inside the same WY%d-%d window sit above their curves by the same margin\n"
                 "unregulated median obs/curve %.2f (%d of %d above)   regulated %.2f (%d of %d above)   "
                 "largest ten unregulated years in window: %d of 10"
                 % (WINDOW + (np.median(u_ratio), int((u_ratio > 1).sum()), n_win,
                              np.median(r_ratio), int((r_ratio > 1).sum()), len(usable), in_win)),
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.95)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print("wrote", OUT_PNG)


if __name__ == "__main__":
    main()
