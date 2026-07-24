#Unreg_Durations_MassBalance.py
# -*- coding: utf-8 -*-
"""
Compute daily unregulated flow at Castle Rock by mass balance and extract
1-, 3-, and 5-day annual maximum durations for the early regulation period.

Method (same as BasicUnregCAS.py, extended to durations):
    Unregulated Flow = Observed Castle Rock Flow + delta(Mossyrock Storage)

Daily mass balance is considered adequate for 1- to 5-day durations but NOT
for instantaneous peak (routing lag and sub-daily storage swings are lost at
daily resolution). Peaks for these years are estimated separately by the
pooled peak-from-1-day regression in Combine_Records.py -- the same
documented procedure already applied to WYs 1982-1987.

Provenance role:
    CDB elevation data for MOS begins in 1974, and the daily Castle Rock
    record on hand covers 11Jan1975 - 28Jan1980. WY1975 and WY1980 are
    flood-season-incomplete (missing Oct-Dec 1974 and Feb-Sep 1980
    respectively), so only WY1976-1979 are adopted from this computation.
    Results are compared against Table B-I of the 2009 Hydrology Restudy
    report ("basic routing model with daily data") as a cross-validation of
    both this computation and the report's method.

Mayfield storage change is not included (INCLUDE_MAY hook below). Mayfield
re-regulation storage is small relative to Mossyrock and largely nets out
over 1-day and longer windows; this is the same simplification made in
BasicUnregCAS.py and should be stated in documentation.
"""

from pydsstools.heclib.dss import HecDss
import pandas as pd
import numpy as np
import os

# Run-from-anywhere: all paths below are relative to this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# =============================================================================
# USER SETTINGS
# =============================================================================
IN_DSS = r"../data/obsData.dss"

MOS_ELEV_PATH    = "//MOS/ELEV-FOREBAY//1DAY/IRVZZAZD_CLEANED/"
CASTLE_ROCK_PATH = "/COWLITZ RIVER AT CASTLE ROCK/14243000/FLOW//1DAY/USGS/"

OUT_CSV = r"../output/unreg_durations_massbalance.csv"

# Water years to adopt (complete flood-season coverage within the data window)
ADOPT_WYS = [1976, 1977, 1978, 1979]

# Report every WY touched by the data, including partials, for transparency
REPORT_ALL_WYS = True

# Mayfield hook -- off to match BasicUnregCAS.py. If enabled, supply a MAY
# elevation path and a MAY elevation-storage table below.
INCLUDE_MAY   = False
MAY_ELEV_PATH = ""  # e.g. "//MAY/ELEV-FOREBAY//1DAY/IRVZZAZD_CLEANED/"

# --- Diagnostic / cleaning ---------------------------------------------------
# Run first with PLOT_DIAGNOSTIC = True and SMOOTH_RANGES empty: inspect the
# four-panel plot and the flagged-days table, identify the erroneous span,
# then add it to SMOOTH_RANGES and re-run. Every smoothing edit is logged to
# the console and must be carried into the memo's data-editing documentation.
PLOT_DIAGNOSTIC = True
DIAG_PLOT_HTML = r"../diagnostics/unreg_massbalance_diagnostic.html"

# Elevation spans to repair by linear interpolation across the range
# (endpoints exclusive of the bad data, i.e. values strictly inside the range
# are replaced). Format: ("ddMMMyyyy", "ddMMMyyyy").
SMOOTH_RANGES = [
    # Non-physical MOS forebay excursion 11-14 Dec 1976 (onset of the 1976-77
    # drought; no water available for a 6-8 ft swing). Without this edit the
    # artifact fabricates a WY1977 1-day max of ~36,200 cfs; with it, WY1977 =
    # 14,719 cfs (Apr 1977), consistent with Table B-I (13,929). VERIFIED
    # 14 Jul 2026: this range reproduces the adopted CSV to 0.000000 cfs for
    # all adopted WYs, and WY1976/1978/1979 are identical with or without it.
    ("01Nov1976", "15Dec1976"),
]

# Flag thresholds for the diagnostic table
ELEV_JUMP_FT_FLAG = 6.0        # |day-over-day forebay change| considered suspect
DELTA_STOR_CFS_FLAG = 60_000   # |storage change| in cfs considered suspect

# Table B-I, 2009 Hydrology Restudy (COWLITZ_HYDROLOGY_REPORT_DRAFT2.docx),
# unregulated 1-day flow (cfs), "basic routing model with daily data".
# !! VERIFY these against the native report table before relying on the
# !! comparison -- values below were read from a rasterized copy.
TABLE_B1_ONEDAY = {
    1976: 114186,
    1977: 13929,
    1978: 115990,
    1979: 33380,
}

# =============================================================================
# CONSTANTS
# =============================================================================
AF_PER_CFS_DAY = 24 * 60 * 60 / 43560.0   # acre-feet per (cfs for one day)

DURATIONS = {          # label -> rolling window in days
    "One_day":   1,
    "Three_Day": 3,
    "Five_Day":  5,
}

FLOOD_SEASON_MONTHS = [10, 11, 12, 1, 2, 3]   # Oct-Mar

# Mossyrock elevation-storage, Table 2-3, 2014 Water Control Manual
# (same table as MOS_CDB_INFLOW.py and BasicUnregCAS.py)
_ELEV = np.array([
    622.0, 623.0, 624.0, 625.0, 626.0, 627.0, 628.0, 629.0, 630.0,
    631.0, 632.0, 633.0, 634.0, 635.0, 636.0, 637.0, 638.0, 639.0, 640.0,
    641.0, 642.0, 643.0, 644.0, 645.0, 646.0, 647.0, 648.0, 649.0, 650.0,
    651.0, 652.0, 653.0, 654.0, 655.0, 656.0, 657.0, 658.0, 659.0, 660.0,
    661.0, 662.0, 663.0, 664.0, 665.0, 666.0, 667.0, 668.0, 669.0, 670.0,
    671.0, 672.0, 673.0, 674.0, 675.0, 676.0, 677.0, 678.0, 679.0, 680.0,
    681.0, 682.0, 683.0, 684.0, 685.0, 686.0, 687.0, 688.0, 689.0, 690.0,
    691.0, 692.0, 693.0, 694.0, 695.0, 696.0, 697.0, 698.0, 699.0, 700.0,
    701.0, 702.0, 703.0, 704.0, 705.0, 706.0, 707.0, 708.0, 709.0, 710.0,
    711.0, 712.0, 713.0, 714.0, 715.0, 716.0, 717.0, 718.0, 719.0, 720.0,
    721.0, 722.0, 723.0, 724.0, 725.0, 726.0, 727.0, 728.0, 729.0, 730.0,
    731.0, 732.0, 733.0, 734.0, 735.0, 736.0, 737.0, 738.0, 739.0, 740.0,
    741.0, 742.0, 743.0, 744.0, 745.0, 746.0, 747.0, 748.0, 749.0, 750.0,
    751.0, 752.0, 753.0, 754.0, 755.0, 756.0, 757.0, 758.0, 759.0, 760.0,
    761.0, 762.0, 763.0, 764.0, 765.0, 766.0, 767.0, 768.0, 769.0, 770.0,
    771.0, 772.0, 773.0, 774.0, 775.0, 776.0, 777.0, 778.0, 779.0, 780.0,
])

_STOR = np.array([
      2_579,   7_736,  12_893,  18_050,  23_405,  28_562,  33_719,  38_876,  44_033,
     49_785,  55_339,  61_091,  66_645,  72_397,  77_950,  83_702,  89_256,  95_008, 100_562,
    106_512, 112_661, 118_612, 124_760, 130_711, 136_661, 142_810, 148_760, 154_711, 160_859,
    167_207, 173_554, 179_901, 186_446, 192_793, 199_140, 205_487, 211_835, 218_380, 224_727,
    231_471, 238_413, 245_157, 251_901, 258_843, 265_587, 272_529, 279_273, 286_016, 292_959,
    300_297, 307_636, 314_777, 322_116, 329_454, 336_793, 344_132, 351_471, 358_810, 366_149,
    373_884, 381_620, 389_355, 397_091, 404_826, 412_562, 420_099, 427_834, 435_570, 443_306,
    451_438, 459_570, 467_702, 475_636, 483_768, 491_901, 500_033, 507_967, 516_099, 524_231,
    532_760, 541_091, 549_620, 558_148, 566_677, 575_008, 583_537, 592_066, 600_595, 608_925,
    617_851, 626_777, 635_504, 644_429, 653_157, 662_082, 671_008, 679_735, 688_661, 697_388,
    706_710, 715_834, 725_157, 734_281, 743_405, 752_727, 761_851, 770_975, 780_297, 789_421,
    798_942, 808_661, 818_181, 827_702, 837_223, 846_743, 856_264, 865_983, 875_504, 885_024,
    894_942, 905_057, 914_975, 924_892, 934_809, 944_925, 954_842, 964_760, 974_677, 984_793,
    995_107, 1_005_619, 1_016_132, 1_026_446, 1_036_958, 1_047_272, 1_057_785, 1_068_297, 1_078_611, 1_089_123,
  1_100_032, 1_111_140, 1_122_049, 1_133_156, 1_144_065, 1_155_173, 1_166_082, 1_177_189, 1_188_099, 1_199_206,
  1_210_710, 1_222_413, 1_233_917, 1_245_619, 1_257_322, 1_268_826, 1_280_528, 1_292_231, 1_303_735, 1_315_437,
], dtype=float)


def read_dss_record(dss_file, pathname, value_name):
    with HecDss.Open(dss_file, version=6) as dss:
        ts = dss.read_ts(pathname)

    df = pd.DataFrame({
        "DateTime": pd.to_datetime(ts.pytimes),
        value_name: ts.values,
    })

    df[value_name] = pd.to_numeric(df[value_name], errors="coerce")
    # -900 threshold: catches -901/-902 DSS sentinels (a -9000 threshold
    # passes them through as data)
    df = df[df[value_name] > -900]
    df = df.dropna()
    df = df.set_index("DateTime")
    return df


def to_water_year(idx):
    return idx.year + (idx.month >= 10).astype(int)


# =============================================================================
# 1. Read records
# =============================================================================
# Guard: pydsstools silently CREATES an empty DSS file if the path does not
# exist, and reads against that empty file fail with a cryptic
# "zcheckKeys ... primary table ... corrupt" error. Fail clearly instead.
if not os.path.isfile(IN_DSS):
    raise FileNotFoundError(
        f"IN_DSS not found: {IN_DSS}\n"
        "Fix the path before running. If a previous run created a small "
        "stray obsData.dss at a wrong location, delete it.")

with HecDss.Open(IN_DSS, version=6) as _dss:
    _catalog = _dss.getPathnameList("/*/*/*/*/*/*/", sort=1)
for _p in (MOS_ELEV_PATH, CASTLE_ROCK_PATH):
    _pat = _p.upper().split("/")
    _hits = [c for c in _catalog
             if all(seg in c.upper() for seg in _pat if seg)]
    if not _hits:
        _preview = "\n  ".join(_catalog[:40])
        raise KeyError(
            f"No record matching {_p} in {IN_DSS}.\n"
            f"First catalog entries:\n  {_preview}")

mos = read_dss_record(IN_DSS, MOS_ELEV_PATH, "MOS_ELEV")
cr  = read_dss_record(IN_DSS, CASTLE_ROCK_PATH, "CASTLE_ROCK_FLOW")

# Normalize to date resolution so joins are exact
mos.index = mos.index.normalize()
cr.index  = cr.index.normalize()

# =============================================================================
# 2. Elevation cleaning (documented edits), then storage change
# =============================================================================
mos["MOS_ELEV_RAW"] = mos["MOS_ELEV"].copy()

for start_s, end_s in SMOOTH_RANGES:
    t0 = pd.to_datetime(start_s, format="%d%b%Y")
    t1 = pd.to_datetime(end_s, format="%d%b%Y")
    inside = (mos.index > t0) & (mos.index < t1)
    n_edit = int(inside.sum())
    if n_edit == 0:
        print(f"EDIT LOG: {start_s}-{end_s}: no interior daily values found "
              f"-- check the range.")
        continue
    if t0 not in mos.index or t1 not in mos.index:
        raise ValueError(
            f"SMOOTH_RANGES endpoints must exist in the elevation record: "
            f"{start_s}, {end_s}")
    before = mos.loc[inside, "MOS_ELEV"].copy()
    span_days = (t1 - t0).days
    frac = (mos.index[inside] - t0).days / span_days
    mos.loc[inside, "MOS_ELEV"] = (
        mos.at[t0, "MOS_ELEV"] + frac * (mos.at[t1, "MOS_ELEV"] - mos.at[t0, "MOS_ELEV"])
    )
    print(f"EDIT LOG: linearly interpolated MOS forebay elevation across "
          f"{start_s} -> {end_s} ({n_edit} interior days replaced; raw range "
          f"{before.min():.2f}-{before.max():.2f} ft, anchors "
          f"{mos.at[t0, 'MOS_ELEV']:.2f} / {mos.at[t1, 'MOS_ELEV']:.2f} ft). "
          f"Document this edit in the memo.")

mos["MOS_STOR_AF"]    = np.interp(mos["MOS_ELEV"], _ELEV, _STOR)
mos["DELTA_STOR_AF"]  = mos["MOS_STOR_AF"].diff()
mos["DELTA_STOR_CFS"] = mos["DELTA_STOR_AF"] / AF_PER_CFS_DAY

if INCLUDE_MAY:
    raise NotImplementedError(
        "Supply a MAY elevation-storage table and mirror the MOS block above."
    )

# =============================================================================
# 3. Daily unregulated flow
# =============================================================================
combined = cr.join(mos["DELTA_STOR_CFS"], how="inner").dropna()
combined["UNREG_CFS"] = combined["CASTLE_ROCK_FLOW"] + combined["DELTA_STOR_CFS"]
combined["WY"] = to_water_year(combined.index)

print(f"Daily unreg record: {combined.index[0].date()} -> "
      f"{combined.index[-1].date()}  ({len(combined)} days)")

# =============================================================================
# 3b. Diagnostics: flag suspect days, plot the record
# =============================================================================
elev_jump = mos["MOS_ELEV"].diff().abs()
flag = pd.DataFrame({
    "MOS_ELEV": mos["MOS_ELEV"],
    "elev_jump_ft": elev_jump,
    "DELTA_STOR_CFS": mos["DELTA_STOR_CFS"],
})
flag = flag[(flag["elev_jump_ft"] > ELEV_JUMP_FT_FLAG) |
            (flag["DELTA_STOR_CFS"].abs() > DELTA_STOR_CFS_FLAG)]
flag = flag.loc[flag.index.intersection(combined.index)]
if len(flag):
    print(f"\nFLAGGED DAYS (|elev jump| > {ELEV_JUMP_FT_FLAG} ft or "
          f"|dS| > {DELTA_STOR_CFS_FLAG:,} cfs) within the unreg record:")
    print(flag.to_string(float_format=lambda v: f"{v:,.1f}"))
    print("If any of these are data artifacts rather than real drawdown/"
          "refill, add the span to SMOOTH_RANGES and re-run.")
else:
    print("\nNo days flagged by elevation-jump / storage-change thresholds.")

if PLOT_DIAGNOSTIC:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.03,
        subplot_titles=("MOS forebay elevation (raw = gray, edited = blue)",
                        "Daily storage change (cfs-equivalent)",
                        "Castle Rock observed flow",
                        "Daily unregulated flow (CR + dStor)"))

    fig.add_trace(go.Scatter(x=mos.index, y=mos["MOS_ELEV_RAW"],
                             name="Elev raw", line=dict(color="lightgray", width=2.2)),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=mos.index, y=mos["MOS_ELEV"],
                             name="Elev edited", line=dict(color="royalblue", width=1)),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=mos.index, y=mos["DELTA_STOR_CFS"],
                             name="dStor (cfs)", line=dict(color="darkorange", width=1)),
                  row=2, col=1)
    fig.add_trace(go.Scatter(x=combined.index, y=combined["CASTLE_ROCK_FLOW"],
                             name="CR observed", line=dict(color="seagreen", width=1)),
                  row=3, col=1)
    fig.add_trace(go.Scatter(x=combined.index, y=combined["UNREG_CFS"],
                             name="Daily unreg", line=dict(color="crimson", width=1)),
                  row=4, col=1)

    if len(flag):
        fig.add_trace(go.Scatter(
            x=flag.index, y=flag["MOS_ELEV"], mode="markers",
            name="flagged day",
            marker=dict(color="magenta", size=8, symbol="x")), row=1, col=1)
        for d in flag.index:
            fig.add_vline(x=d, line_width=1, line_color="magenta",
                          opacity=0.35)

    fig.update_layout(
        height=1000, hovermode="x unified",
        title=("Mass-balance diagnostic -- zoom to the flagged spike, note "
               "the anchor dates, add the span to SMOOTH_RANGES and re-run"))
    fig.update_yaxes(title_text="ft", row=1, col=1)
    fig.update_yaxes(title_text="cfs", row=2, col=1)
    fig.update_yaxes(title_text="cfs", row=3, col=1)
    fig.update_yaxes(title_text="cfs", row=4, col=1)

    fig.write_html(DIAG_PLOT_HTML)
    print(f"\nWrote diagnostic plot: {DIAG_PLOT_HTML}")

# =============================================================================
# 4. Per-WY completeness
# =============================================================================
rows = []
for wy, grp in combined.groupby("WY"):
    wy_start = pd.Timestamp(year=wy - 1, month=10, day=1)
    wy_end   = pd.Timestamp(year=wy, month=9, day=30)
    full_days = (wy_end - wy_start).days + 1

    season_expected = pd.date_range(wy_start, wy_end, freq="D")
    season_expected = season_expected[season_expected.month.isin(FLOOD_SEASON_MONTHS)]
    season_present  = grp.index[grp.index.month.isin(FLOOD_SEASON_MONTHS)]

    row = {
        "WY": wy,
        "days_present": len(grp),
        "days_missing": full_days - len(grp),
        "flood_season_missing": len(season_expected) - len(season_present),
    }

    for label, window in DURATIONS.items():
        roll = grp["UNREG_CFS"].rolling(window, min_periods=window).mean()
        row[label] = roll.max()
        row[f"{label}_Date"] = roll.idxmax().date() if roll.notna().any() else None

    row["TableB1_OneDay"] = TABLE_B1_ONEDAY.get(wy, np.nan)
    row["Adopt"] = wy in ADOPT_WYS
    rows.append(row)

result = pd.DataFrame(rows).set_index("WY").sort_index()
result["OneDay_vs_TableB1_pct"] = (
    (result["One_day"] - result["TableB1_OneDay"]) / result["TableB1_OneDay"] * 100
)
result["Source"] = "MassBalance_Daily"

if not REPORT_ALL_WYS:
    result = result.loc[result["Adopt"]]

# =============================================================================
# 5. Report and write
# =============================================================================
cols = ["days_present", "days_missing", "flood_season_missing",
        "One_day", "One_day_Date", "Three_Day", "Three_Day_Date",
        "Five_Day", "Five_Day_Date",
        "TableB1_OneDay", "OneDay_vs_TableB1_pct", "Adopt"]

print("\nWY duration maxima (cfs):")
print(result[cols].to_string(float_format=lambda v: f"{v:,.0f}"))

partial = result[(result["Adopt"]) & (result["flood_season_missing"] > 0)]
if len(partial):
    print("\nWARNING: adopted WYs with missing flood-season days:")
    print(partial[["days_missing", "flood_season_missing"]].to_string())

print("\nNon-adopted WYs are reported for transparency only; WY1975 and "
      "WY1980 are flood-season-incomplete and must not enter the record.")

result.to_csv(OUT_CSV)
print(f"\nWrote: {OUT_CSV}")
print("\nNext step: peak estimates for adopted WYs come from the pooled "
      "peak-from-1-day regression (Combine_Records.py), matching the "
      "procedure used for WYs 1982-1987.")
