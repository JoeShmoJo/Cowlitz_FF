#Transform_Convergence_Figure.py
# -*- coding: utf-8 -*-
"""
Regulated against unregulated peak at Castle Rock, with the transform drawn
converging onto the 1:1 line.

THIS FIGURE IS A DRAWING, NOT A RESULT
    It writes a PNG and nothing else. No CSV, no curve ordinates, nothing any
    other script reads. The convergence limb is a hand-set estimate and must
    not be pulled into the frequency curves, the downstream locations, or the
    memo tables. Everything left of the last frequency-curve ordinate is the
    adopted transform and is real; everything right of it is illustration.

WHY THE TRANSFORM CONVERGES ON 1:1 AT ALL
    A reservoir cannot remove volume from a flood, only move it in time (see
    the memo's Section 5.2 -- at the 5-day duration the regulated/unregulated
    ratio already exceeds one). Riffe Lake holds 358,116 acre-feet between
    the winter rule curve at 745.5 ft and full pool at 778.5 ft, which is
    180,550 cfs-days. Once an event's pre-crest inflow fills that, every
    further cfs passes through and the regulated peak equals the unregulated
    one.

    The adopted transform already shows this turning. Its maximum reduction
    is 67,558 cfs near the 0.2% event and it DECLINES beyond that:

        AEP      unreg      reg     reduction
        0.005  195,977  128,798      67,179
        0.002  227,900  160,342      67,558   <- maximum
        0.001  255,074  189,434      65,640
        0.0002 332,739  282,959      49,781
        0.0001 374,643  339,018      35,625

    That maximum is a useful check on itself: 180,550 cfs-days divided by
    67,558 cfs is 2.7 days, a plausible time from flood onset to crest in
    this basin. The reduction falls off past that point because larger events
    fill the pool sooner, which is exactly the mechanism that ends in 1:1.

WHY THE CONVERGENCE POINT IS DRAWN RATHER THAN FITTED
    It was fitted, and the fit does not survive contact with the data. The
    48 synthetic members cluster between 265,000 and 280,000 cfs unregulated
    and their reductions there span 1,870 to 99,478 cfs -- Dec1933 passes 99
    percent of its inflow while Dec2015, at almost the same magnitude, is
    still holding back 98,000 cfs. Convergence depends on hydrograph shape
    and starting pool, not on peak magnitude alone, so a regression through
    those points is meaningless: fitting all 48 puts the crossing at 931,600
    cfs, fitting the top twelve puts it at 259,600. Nothing in between is
    better evidenced.

    CONVERGE_AT_CFS below is therefore a judgment, set to be consistent with
    the adopted transform's own declining reduction, which extrapolates to
    zero near 480,000-500,000 cfs. It is one constant, it is labelled on the
    figure as an estimate, and it should be moved if a reviewer prefers a
    different value. Scaling more events would settle it properly and is the
    honest alternative if this ever needs to be defended as a number.

OUTPUT
    ../output/diagnostics/transform_convergence.png     (and nothing else)
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
REG_CSV = r"../output/regulated_frequency_inferred.csv"
HIST_CSV = r"../output/adjusted_peaks.csv"
SYNTH_CSV = r"../output/diagnostics/ResSim_Synth_reg_vs_unreg_wy.csv"
PLOT_PNG = r"../output/diagnostics/transform_convergence.png"

CONVERGE_AT_CFS = 350000.0    # THE DRAWN ESTIMATE. See the docstring.
# Draw the limb from the last SUPPORTED ordinate rather than the last ordinate
# of all. The adopted transform is clipped at the 1:1 line, so its final
# extrapolated ordinate already sits on 1:1 and a limb starting there would
# have nothing to draw. Starting from the last supported point makes the drawn
# estimate a genuine alternative to the curve's own extrapolation, which is
# what the reader needs to see.
LIMB_FROM_LAST_SUPPORTED = True
TAPER_POWER = 1.6             # shape of the drawn limb only; 1 is a straight
                              # run-in, higher lands on 1:1 more gently.

FLOOD_STORAGE_ACFT = 358116.0   # Riffe, 745.5 -> 778.5 ft
AXIS_MAX_CFS = 430000.0

C_1TO1 = "#8a8a8a"
# DQC: Section 5.4's two figures are merged into this one, so the event pairs
# that place the transform are drawn again. Set False for a lines-only figure.
SHOW_EVENT_POINTS = True
# Transform uncertainty, drawn from the sigma the frequency script wrote.
SHOW_TRANSFORM_BAND = True
TRANSFORM_BAND_Z = 1.645        # 90% two-sided (5%/95%), matching Section 5.6

C_CURVE = "#1a4f8a"
C_DRAWN = "#b7410e"
C_HIST = "#4c8c4a"
C_SYNTH = "#d99b30"

# ----------------------------------------------------------------------------


def drawn_limb(u_last, reg_last, converge_at, n=160):
    """Smooth run from the last real ordinate onto the 1:1 line.

    The reduction is tapered to zero at converge_at; the curve is therefore
    monotone and lands tangentially rather than kinking into 1:1.
    """
    shave_last = u_last - reg_last
    u = np.linspace(u_last, converge_at, n)
    frac = (u - u_last) / (converge_at - u_last)
    shave = shave_last * (1.0 - frac) ** TAPER_POWER
    return u, u - shave


def main():
    reg = pd.read_csv(REG_CSV).sort_values("unreg_expected_cfs")
    is_ex = reg["extrapolated"].astype(bool)
    real = reg[~is_ex]
    # Include the last supported ordinate so the two segments share a point
    # and the line does not break where the styling changes.
    extrap = reg[is_ex.shift(-1, fill_value=False) | is_ex]

    # The event pairs are NOT plotted -- this figure is lines only. They are
    # still read so the run prints the scatter that makes the convergence
    # point a judgment rather than a fit; see the docstring.
    hist = pd.read_csv(HIST_CSV)
    if "screen_passed" in hist.columns:
        hist = hist[hist["screen_passed"].astype(bool)]
    hist = hist[["unreg_ref", "adjusted_peak"]].dropna()
    synth = pd.read_csv(SYNTH_CSV)[["unreg_peak", "reg_peak"]].dropna()

    if LIMB_FROM_LAST_SUPPORTED and len(real):
        u_last = float(real["unreg_expected_cfs"].iloc[-1])
        reg_last = float(real["reg_inferred_cfs"].iloc[-1])
    else:
        u_last = float(reg["unreg_expected_cfs"].iloc[-1])
        reg_last = float(reg["reg_inferred_cfs"].iloc[-1])
    if CONVERGE_AT_CFS <= u_last:
        raise SystemExit(
            "CONVERGE_AT_CFS (%s) is at or below the last supported ordinate "
            "(%s).\n  The limb would run backwards. Raise it."
            % (format(int(CONVERGE_AT_CFS), ","), format(int(u_last), ",")))
    u_draw, reg_draw = drawn_limb(u_last, reg_last, CONVERGE_AT_CFS)

    fig, ax = plt.subplots(figsize=(9.5, 8.6))
    lim = [0, AXIS_MAX_CFS]
    ax.plot(lim, lim, color=C_1TO1, lw=1.6, ls="--", zorder=1,
            label="1:1  (no regulation effect)")

    # The transform's own uncertainty, which Section 5.6 combines with the
    # frequency sigma. Drawn here so the reader can see where it came from.
    if SHOW_TRANSFORM_BAND and {"sigma_transform_lo_dex",
                                "sigma_transform_hi_dex"} <= set(reg.columns):
        base = reg["reg_inferred_cfs"].values
        lo = base / 10.0 ** (TRANSFORM_BAND_Z
                             * reg["sigma_transform_lo_dex"].values)
        hi = base * 10.0 ** (TRANSFORM_BAND_Z
                             * reg["sigma_transform_hi_dex"].values)
        ax.fill_between(reg["unreg_expected_cfs"].values, lo, hi,
                        color=C_CURVE, alpha=0.13, lw=0, zorder=2,
                        label="Transform scatter, 5%-95%")

    if SHOW_EVENT_POINTS:
        ax.plot(hist["unreg_ref"], hist["adjusted_peak"], ls="none",
                marker="o", ms=4.6, mfc="none", mec=C_HIST, mew=1.0,
                zorder=3, label="Observed events, adjusted (n=%d)" % len(hist))
        ax.plot(synth["unreg_peak"], synth["reg_peak"], ls="none",
                marker="s", ms=4.6, mfc="none", mec=C_SYNTH, mew=1.0,
                zorder=3, label="Synthetic members (n=%d)" % len(synth))

    ax.plot(real["unreg_expected_cfs"], real["reg_inferred_cfs"],
            color=C_CURVE, lw=2.8, zorder=4, label="Adopted transform")
    ax.plot(extrap["unreg_expected_cfs"], extrap["reg_inferred_cfs"],
            color=C_CURVE, lw=2.8, ls=(0, (6, 2)), zorder=4,
            label="Adopted transform, extrapolated")
    ax.plot(u_draw, reg_draw, color=C_DRAWN, lw=2.4, ls=(0, (2, 2)), zorder=4,
            label="Drawn convergence, straight from the last supported point")

    ax.plot([CONVERGE_AT_CFS], [CONVERGE_AT_CFS], marker="o", ms=11,
            mfc="none", mec=C_DRAWN, mew=2.4, zorder=5)
    ax.annotate("convergence\n~%s cfs" % format(int(CONVERGE_AT_CFS), ","),
                xy=(CONVERGE_AT_CFS, CONVERGE_AT_CFS),
                xytext=(CONVERGE_AT_CFS + 18000, CONVERGE_AT_CFS - 88000),
                color=C_DRAWN, fontsize=9.5, ha="left",
                arrowprops=dict(arrowstyle="->", color=C_DRAWN, lw=1.3))

    # note = ("Riffe usable flood storage %s ac-ft = %s cfs-days\n"
    #         "Maximum reduction on the adopted transform %s cfs,\n"
    #         "which fills that storage in %.1f days."
    #         % (format(int(FLOOD_STORAGE_ACFT), ","),
    #            format(int(FLOOD_STORAGE_ACFT / 1.98347), ","),
    #            format(int((reg["unreg_expected_cfs"]
    #                        - reg["reg_inferred_cfs"]).max()), ","),
    #            (FLOOD_STORAGE_ACFT / 1.98347)
    #            / (reg["unreg_expected_cfs"] - reg["reg_inferred_cfs"]).max()))
    # ax.text(0.035, 0.965, note, transform=ax.transAxes, va="top", fontsize=8.5,
    #         bbox=dict(boxstyle="round,pad=0.5", fc="#f4f1ea", ec="#c9c2b4"))

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Unregulated peak at Castle Rock (cfs)")
    ax.set_ylabel("Regulated peak at Castle Rock (cfs)")
    ax.set_title("Regulated against unregulated peak at Castle Rock\n"
             "adopted transform, its scatter, and the drawn convergence "
             "on 1:1", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.0, framealpha=0.92)
    ax.xaxis.set_major_formatter(lambda v, _: format(int(v / 1000), ",") + "k")
    ax.yaxis.set_major_formatter(lambda v, _: format(int(v / 1000), ",") + "k")

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("Wrote %s" % PLOT_PNG)
    print("   drawn convergence at %s cfs -- a judgment, one constant, "
          "used nowhere else" % format(int(CONVERGE_AT_CFS), ","))
    top = synth[synth["unreg_peak"] >= 260000]
    shave = top["unreg_peak"] - top["reg_peak"]
    print("   lines only; %d historical and %d synthetic pairs read but not "
          "plotted" % (len(hist), len(synth)))
    print("   why it is drawn: %d synthetic members above 260,000 cfs "
          "unregulated\n   span reductions of %s to %s cfs -- shape and "
          "starting pool, not magnitude"
          % (len(top), format(int(shave.min()), ","),
             format(int(shave.max()), ",")))


if __name__ == "__main__":
    main()
