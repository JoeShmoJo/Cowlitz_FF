#Coweeman_PRISM_PrecipRatio.py
# -*- coding: utf-8 -*-
"""
Coweeman / Cowlitz-above-Castle-Rock basin precipitation ratio from PRISM.

WHY
    The peak ratio in #BelowConfluence_FlowFrequency.py is a drainage-area
    ratio: it assumes the two basins receive and convert precipitation at the
    same rate per square mile. This tests that assumption directly. If the
    Coweeman's basin-mean precipitation is close to Castle Rock's, an area
    ratio is defensible on its own terms; if it runs high or low, that offset
    is the physical explanation for the flow ratio sitting off 1.00x, and it
    transfers to the ungaged Arkansas and Ostrander basins the same way.

    Annual totals are enough. The Cowlitz's annual maximum is a winter
    frontal/rain-on-snow event and those storms dominate the annual total, so
    the annual ratio tracks the storm ratio closely enough for this purpose.
    This is a check on the AREA assumption, not a rainfall-runoff model.

WHY THIS SCRIPT DOES NOT RUN IN THE SANDBOX
    It needs two things that must come from outside: PRISM grids
    (services.nacse.org) and basin polygons. Both are blocked from the remote
    session -- every data host tested returned a connection failure while
    pypi.org returned 200, so it is host blocking, not an outage. Run this on
    your own machine.

WHAT YOU MUST SUPPLY
    Two basin polygons, in BASIN_FILES below. Any format geopandas reads
    (.shp, .geojson, .gpkg). Sources, easiest first:

      1. USGS StreamStats (streamstats.usgs.gov) -- delineate from the gage
         location and export the basin. Coweeman at USGS 14245000; Cowlitz at
         USGS 14243000 (Castle Rock).
      2. NHDPlus HR catchment aggregation, or the WBD HUC-10 for the Coweeman
         (1708000504) if a gage-exact boundary is not needed.

    Check the delineated areas against the numbers this study uses --
    Coweeman 119 sq mi at the gage, Cowlitz 2,238 sq mi above Castle Rock.
    The script prints polygon areas so a bad boundary shows up immediately
    rather than silently biasing the ratio.

WHAT IT DOES
    Downloads PRISM annual precipitation (4km, ppt, stable), clips each grid
    to each basin, takes the area-weighted mean, and reports the ratio by
    year. all_touched=False so a cell counts only where its centre falls
    inside the basin -- on a 4km grid against a 119 sq mi basin that matters:
    all_touched=True inflates a small basin by grabbing partial edge cells.

    A 119 sq mi basin is only about 7-8 PRISM cells. The Coweeman mean is
    therefore coarse, and the script prints the cell count so that is visible.
    If it is under MIN_CELLS_WARN the ratio should be treated as indicative.

OUTPUTS
    ../output/diagnostics/prism_basin_precip_ratio.csv
    ../output/diagnostics/prism_basin_precip_ratio.png
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import glob
import io
import zipfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
BASIN_FILES = {
    "Coweeman": r"../data/basins/coweeman_14245000.geojson",
    "CastleRock": r"../data/basins/cowlitz_14243000.geojson",
}
EXPECTED_SQ_MI = {"Coweeman": 119.0, "CastleRock": 2238.0}

PRISM_DIR = r"../data/prism"          # grids cached here
YEAR_START, YEAR_END = 1950, 2020     # PRISM stable series starts 1895
PRISM_URL = "https://services.nacse.org/prism/data/public/4km/ppt/%d"

OUT_CSV = r"../output/diagnostics/prism_basin_precip_ratio.csv"
PLOT_PNG = r"../output/diagnostics/prism_basin_precip_ratio.png"

MIN_CELLS_WARN = 10        # below this the basin mean is coarse; say so
NODATA = -9999.0

C_COW = "#b7410e"
C_CAS = "#1a4f8a"

# ----------------------------------------------------------------------------


def need(module):
    try:
        return __import__(module)
    except ImportError:
        raise SystemExit(
            "%s is required.  pip install geopandas rasterio requests" % module)


def fetch_year(year):
    """Download and unpack one PRISM annual ppt grid, cached."""
    import requests
    os.makedirs(PRISM_DIR, exist_ok=True)
    hit = glob.glob(os.path.join(PRISM_DIR, "*_%d_*.bil" % year))
    if hit:
        return hit[0]
    url = PRISM_URL % year
    print("   downloading %d ..." % year, end="", flush=True)
    resp = requests.get(url, timeout=180)
    resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as archive:
        archive.extractall(PRISM_DIR)
    hit = glob.glob(os.path.join(PRISM_DIR, "*_%d_*.bil" % year))
    if not hit:
        raise SystemExit("no .bil found in the %d archive" % year)
    print(" ok")
    return hit[0]


def load_basins():
    gpd = need("geopandas")
    out = {}
    for name, path in BASIN_FILES.items():
        if not os.path.exists(path):
            raise SystemExit(
                "Basin polygon not found: %s\nSee 'WHAT YOU MUST SUPPLY' in "
                "this script's docstring." % os.path.abspath(path))
        frame = gpd.read_file(path)
        # Albers equal area for an honest area check; PRISM itself is
        # geographic, so the clip happens in the raster's own CRS later.
        area = frame.to_crs("EPSG:5070").area.sum() / 2.589988e6
        want = EXPECTED_SQ_MI[name]
        flag = "" if abs(area - want) / want < 0.10 else "   <-- CHECK BOUNDARY"
        print("   %-12s %8.1f sq mi  (study uses %.0f)%s" % (name, area, want, flag))
        out[name] = frame
    return out


def basin_mean(raster_path, basin):
    """Area-weighted mean of a PRISM grid over one polygon, plus cell count."""
    rasterio = need("rasterio")
    from rasterio.mask import mask
    with rasterio.open(raster_path) as src:
        shapes = basin.to_crs(src.crs).geometry.values
        clipped, _ = mask(src, shapes, crop=True, filled=True,
                          nodata=NODATA, all_touched=False)
    values = clipped[0].astype(float)
    values[values == NODATA] = np.nan
    good = np.isfinite(values)
    if not good.any():
        return np.nan, 0
    return float(np.nanmean(values)), int(good.sum())


def main():
    print("Basins:")
    basins = load_basins()

    rows = []
    print("\nPRISM annual precipitation:")
    for year in range(YEAR_START, YEAR_END + 1):
        try:
            grid = fetch_year(year)
        except Exception as exc:                     # noqa: BLE001
            print("   %d skipped: %s" % (year, exc))
            continue
        row = {"year": year}
        for name, basin in basins.items():
            mean_mm, cells = basin_mean(grid, basin)
            row["%s_ppt_mm" % name] = mean_mm
            row["%s_cells" % name] = cells
        row["ratio"] = (row["Coweeman_ppt_mm"] / row["CastleRock_ppt_mm"]
                        if row["CastleRock_ppt_mm"] else np.nan)
        rows.append(row)

    table = pd.DataFrame(rows)
    if not len(table):
        raise SystemExit("No PRISM years were retrieved.")
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    table.to_csv(OUT_CSV, index=False)

    cells = int(table["Coweeman_cells"].median())
    print("\n" + "=" * 74)
    print("BASIN-MEAN ANNUAL PRECIPITATION, %d-%d" % (table.year.min(), table.year.max()))
    print("=" * 74)
    print("   Coweeman   %7.0f mm   (%d PRISM cells%s)"
          % (table["Coweeman_ppt_mm"].mean(), cells,
             ", COARSE" if cells < MIN_CELLS_WARN else ""))
    print("   CastleRock %7.0f mm   (%d PRISM cells)"
          % (table["CastleRock_ppt_mm"].mean(),
             int(table["CastleRock_cells"].median())))
    r = table["ratio"]
    print("\n   precipitation ratio Coweeman / CastleRock")
    print("      mean %.3f   median %.3f   p25 %.3f   p75 %.3f   min %.3f   max %.3f"
          % (r.mean(), r.median(), r.quantile(.25), r.quantile(.75), r.min(), r.max()))

    print("\n   HOW TO READ IT")
    print("      ratio ~ 1.00  the two basins get the same depth, so a plain")
    print("                    drainage-area ratio is defensible.")
    print("      ratio > 1.00  the Coweeman is wetter per unit area, which")
    print("                    would explain a flow ratio above 1.00x area.")
    print("      ratio < 1.00  drier, and an area ratio OVERSTATES the local.")
    print("\n      Compare against the flow evidence in")
    print("      #Coweeman_HistoricPeakRatio.py: the Coweeman peak runs about")
    print("      1.5x its area share at common events but converges to about")
    print("      1.0x at the largest ones. Precipitation that is roughly equal")
    print("      would say the common-event excess is a RESPONSE difference")
    print("      (small flashy basin, faster concentration) rather than a")
    print("      precipitation difference -- and response differences shrink")
    print("      as both basins saturate, which is what the flow data shows.")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(table.year, table["Coweeman_ppt_mm"], color=C_COW, lw=1.6, label="Coweeman")
    ax1.plot(table.year, table["CastleRock_ppt_mm"], color=C_CAS, lw=1.6, label="Castle Rock")
    ax1.set_xlabel("year")
    ax1.set_ylabel("basin-mean annual precipitation (mm)")
    ax1.set_title("PRISM basin-mean annual precipitation")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=9)

    ax2.plot(table.year, table["ratio"], color="#4c8c4a", lw=1.6)
    ax2.axhline(1.0, color="0.4", ls="--", lw=1.4, label="equal depth")
    ax2.axhline(r.mean(), color="#b7410e", ls=":", lw=1.8,
                label="mean %.3f" % r.mean())
    ax2.set_xlabel("year")
    ax2.set_ylabel("Coweeman / Castle Rock precipitation")
    ax2.set_title("Precipitation ratio\nthe physical test of the area-ratio assumption")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print("\nWrote", OUT_CSV)
    print("Wrote", PLOT_PNG)


if __name__ == "__main__":
    main()
