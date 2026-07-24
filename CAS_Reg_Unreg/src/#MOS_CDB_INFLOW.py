#MOS_CDB_INFLOW.py
'''
Calculate MOS inflow from hourly elevation and outflow using mass balance.
Uses MOS elevation - storage table from 2014 Water Control Manual
Needed for old data, such as the CDB records, which have MOS elevation but not storage.
Used here to recreate the MOS inflow time series for 1996, which is needed to test 
Emergency Spillway Releases in ResSim.
'''

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Settings — edit these
# ---------------------------------------------------------------------------
CSV_IN  = "../data/MOS_96.csv"
CSV_OUT = "../output/MOS_96_recalc.csv"

# ---------------------------------------------------------------------------
# Rating curve: Mossyrock Table 2-3 (Sheet 1 of 2)
# Elev (ft) → Volume (acre-ft)
# Note: 621.0 ft has no listed volume; table starts at 622.0 ft
# ---------------------------------------------------------------------------
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

_interp = PchipInterpolator(_ELEV, _STOR)

def elev_to_stor(elev):
    return _interp(elev)

# ---------------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------------
df = pd.read_csv(CSV_IN, skipinitialspace=True)
df.columns = df.columns.str.strip()

# Handle trailing comma on Year field e.g. "1996,"
df["Year"] = df["Year"].astype(str).str.replace(r"[,\s]+$", "", regex=True).astype(int)

df["datetime"] = pd.to_datetime(
    df["Year"].astype(str) + "-" + df["Month"].astype(str) + "-"
    + df["Day"].astype(str) + " " + df["Time"].astype(str),
    format="%Y-%b-%d %H:%M",
)
df = df.sort_values("datetime").reset_index(drop=True)

# Coerce flow columns to numeric (handles comma-formatted strings e.g. "8,200")
for col in ["FLOW-IN", "FLOW-OUT"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col].astype(str).str.replace(",", ""), errors="coerce")

# ---------------------------------------------------------------------------
# Interpolate instantaneous observations to hourly
# ---------------------------------------------------------------------------
df = df.set_index("datetime")

hourly_index = pd.date_range(df.index[0], df.index[-1], freq="h")
df_hr = df[["ELEV", "FLOW-OUT", "FLOW-IN"]].reindex(
    df.index.union(hourly_index)
).interpolate(method="time").reindex(hourly_index)

# ---------------------------------------------------------------------------
# Mass balance at hourly time step
# Inflow = dS/dt + Outflow   (Δt = 1 hr)
# ---------------------------------------------------------------------------
CFS_PER_ACFT_HR = 43_560 / 3_600   # 12.1008

df_hr["STOR_AF"]        = elev_to_stor(df_hr["ELEV"].values)
df_hr["DS_AF"]          = df_hr["STOR_AF"].diff()
df_hr["INFLOW_RECALC"]  = df_hr["DS_AF"] * CFS_PER_ACFT_HR + df_hr["FLOW-OUT"]

df_hr.index.name = "datetime"

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
df_hr.to_csv(CSV_OUT, float_format="%.3f")
print(f"Done. {len(df_hr)} hourly rows written to {CSV_OUT}")
print(df_hr[["ELEV", "STOR_AF", "DS_AF", "FLOW-OUT", "FLOW-IN", "INFLOW_RECALC"]].to_string())

# create MOS stor elevation table at every foot  between 621 and 780 ft
elevs = np.arange(621.0, 780.0 + 1.0, 1.0)
stor_af = elev_to_stor(elevs)
mos_table = pd.DataFrame({"ELEV": elevs, "STOR_AF": stor_af})
print("\nMOS Storage Table:")