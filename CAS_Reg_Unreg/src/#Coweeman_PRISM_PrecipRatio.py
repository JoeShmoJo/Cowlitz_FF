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

GETTING THE GRIDS -- OFFLINE FIRST, ON PURPOSE
    The download is a convenience, not a dependency. PRISM has restructured
    its download service more than once, and services.nacse.org fails to
    resolve on at least one network this has been run from (DNS, errno 11002)
    -- so the script uses whatever grids are already in PRISM_DIR, tries to
    fetch only the years that are missing, and gives up on the network after
    the FIRST failure rather than repeating the same error once per year.

    If the download does not work, download annual 4km precipitation by hand
    from https://prism.oregonstate.edu/ and drop the .zip files straight into
    PRISM_DIR -- the script unpacks them. Filenames do not matter as long as
    the year appears in them; .bil and .tif are both read.

    Partial coverage is fine. The ratio is a stable basin property, not a peak
    statistic, so ten scattered years settle it about as well as seventy.

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
    # StreamStats "globalwatershed" delineations, uploaded and verified:
    # 118.2 and 2229.3 sq mi against StreamStats' own report, within 0.1%.
    "Coweeman": r"../data/BasinDelinations/Coweeman_SHP/layers/globalwatershed.shp",
    "CastleRock": r"../data/BasinDelinations/CastleRock_SHP/layers/globalwatershed.shp",
}
EXPECTED_SQ_MI = {"Coweeman": 118.2, "CastleRock": 2229.3}

PRISM_DIR = r"../data/prism"          # grids live/cache here
YEAR_START, YEAR_END = 1950, 2020     # PRISM stable series starts 1895

# Downloading is OPTIONAL and is tried only for years not already present.
# PRISM has restructured its download service more than once and this host
# does not resolve on every network -- if it fails, the script says exactly
# what to download by hand and carries on with whatever is already local.
# Paste a working URL pattern here if you have one; %d is the year.
# Per CAS_Reg_Unreg/data/PRISM_downloads_web_service.pdf (26 Mar 2025):
# .../data/get/<region>/<res>/<element>/<date>?format=bil -- one grid per
# request, delivered as COG (.tif) by default; ?format=bil asks for the BIL
# variant instead so the rest of this script does not need to change. The
# OLD path this used, .../data/public/4km/ppt/<year>, is retired.
PRISM_URL = "https://services.nacse.org/prism/data/get/us/4km/ppt/%d?format=bil"
ALLOW_DOWNLOAD = True
# The service blocks a second download of the same file within 24h (see the
# PDF's DOWNLOAD LIMITS section) -- another reason gather_years() never
# re-requests a year that is already on disk.

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


def unpack_local_zips():
    """Unpack any PRISM .zip the user dropped in PRISM_DIR by hand."""
    for archive_path in sorted(glob.glob(os.path.join(PRISM_DIR, "*.zip"))):
        try:
            with zipfile.ZipFile(archive_path) as archive:
                archive.extractall(PRISM_DIR)
        except zipfile.BadZipFile:
            print("   not a zip, skipped: %s" % os.path.basename(archive_path))


def local_grid(year):
    """Any already-present grid for this year, .bil or .tif, or None."""
    for pattern in ("*%d*.bil" % year, "*%d*.tif" % year, "*%d*.tiff" % year):
        hit = sorted(glob.glob(os.path.join(PRISM_DIR, pattern)))
        if hit:
            return hit[0]
    return None


def try_download(year):
    """Fetch one year. Returns (path, None) or (None, short reason)."""
    import requests
    try:
        resp = requests.get(PRISM_URL % year, timeout=180)
        resp.raise_for_status()
    except Exception as exc:                                  # noqa: BLE001
        text = str(exc)
        if "getaddrinfo" in text or "NameResolution" in text:
            return None, "DNS: %s does not resolve on this network" % (
                PRISM_URL.split("/")[2])
        return None, text.split("(Caused by")[0].strip()[:120]
    try:
        with zipfile.ZipFile(io.BytesIO(resp.content)) as archive:
            archive.extractall(PRISM_DIR)
    except zipfile.BadZipFile:
        return None, "response was not a zip (URL pattern probably stale)"
    path = local_grid(year)
    return (path, None) if path else (None, "no grid found inside the archive")


def gather_years():
    """Grids for every year we can get, downloading only what is missing.

    Offline-first on purpose. A locked-down network is the normal case here,
    and one clear instruction beats seventy identical stack traces.
    """
    os.makedirs(PRISM_DIR, exist_ok=True)
    unpack_local_zips()
    found, missing, reason = {}, [], None
    for year in range(YEAR_START, YEAR_END + 1):
        path = local_grid(year)
        if path:
            found[year] = path
            continue
        if not ALLOW_DOWNLOAD or reason is not None:
            missing.append(year)          # already know the network is out;
            continue                      # do not retry it 70 more times
        path, reason = try_download(year)
        if path:
            found[year] = path
            reason = None
        else:
            missing.append(year)

    print("   %d year(s) available locally" % len(found))
    if missing:
        print("   %d year(s) missing" % len(missing))
        if reason:
            print("   download unavailable -- %s" % reason)
        print("""
   TO SUPPLY THEM BY HAND
     1. Open  https://prism.oregonstate.edu/  ->  Data Explorer / downloads
        and take ANNUAL total precipitation, 4km, for the years you want.
     2. Drop the .zip files (or the unpacked .bil/.hdr/.prj sets) into
          %s
        Filenames do not matter as long as the year appears in them.
     3. Re-run. Anything already there is used and never re-downloaded.

   Partial coverage is fine -- the ratio is computed per year and averaged
   over whatever is present. Even ten scattered years settles this question,
   since the ratio is a stable basin property rather than a peak statistic.
""" % os.path.abspath(PRISM_DIR))
    return found


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

    print("\nPRISM annual precipitation:")
    grids = gather_years()
    if not grids:
        raise SystemExit(
            "No PRISM grids available. See the instructions above -- the "
            "basin polygons loaded fine, so this is the only thing missing.")

    rows = []
    for year, grid in sorted(grids.items()):
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
