#make_basin_map.py
# -*- coding: utf-8 -*-
"""
Cowlitz River basin location map, on an ESRI basemap.

Replaces basin_schematic.png, which was a hand-drawn matplotlib cartoon with
made-up coordinates and "not to scale" in the title. Everything here is real
geography: the river lines are NHDPlus flowlines from the USGS Network Linked
Data Index (NLDI), the gage positions are the coordinates NWIS publishes for
those station numbers, and the backdrop is an ESRI tile service.

WHAT IS DRAWN
    - the Cowlitz mainstem, from Castle Rock up to the headwaters
    - the Toutle, drawn as its own line because it is the tributary that
      matters here: it enters BELOW both dams, so its flow is unregulated and
      lands on the Castle Rock gage unattenuated
    - the rest of the upstream network, thin, for context
    - the drainage basin above Castle Rock
    - Mossyrock and Mayfield dams, and Riffe / Mayfield lakes behind them
    - the gages: Castle Rock (14243000) and the Mayfield outflow (14238000)

WHY IT NEEDS ITS OWN ENVIRONMENT
    geopandas / contextily / rasterio pull a modern numpy, and pyogrio and
    shapely 2 are strict about it. Dropping them into the analysis environment
    risks moving numpy and pandas underneath the DSS and frequency scripts,
    which is not a trade worth making for a figure that is generated once.
    environment.yml beside this file builds a standalone env; nothing in the
    analysis chain imports this script.

        conda env create -f environment.yml
        conda activate cowlitz-map
        python make_basin_map.py

NETWORK, AND THE CACHE
    The first run fetches from NLDI and NWIS and writes what it gets to
    CACHE_DIR as GeoJSON/CSV. Later runs read the cache and need no network
    except for the basemap tiles, which contextily caches separately. Delete
    CACHE_DIR to force a refresh.

    If a fetch fails the script says which one and carries on with whatever it
    has, so a network problem degrades the map rather than killing it. It will
    NOT silently invent geometry: FALLBACK_SITES exists only so the gage
    markers can still be placed, and every fallback is announced and stamped on
    the figure.
"""

import os
# Run-from-anywhere: relative paths below resolve from this script's folder
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import json
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import contextily as cx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from shapely.geometry import Point

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
OUT_PNG = r"basin_map.png"
CACHE_DIR = r"mapdata"
DPI = 300
FIGSIZE = (12.0, 8.0)

# --- basemap -----------------------------------------------------------------
# Any of: WorldTopoMap, WorldImagery, WorldTerrain, WorldShadedRelief,
# WorldStreetMap, WorldPhysical, NatGeoWorldMap, WorldGrayCanvas.
#   WorldTopoMap   labelled topo -- place names and roads. The default: it
#                  already draws and labels Riffe and Mayfield lakes.
#   WorldImagery   aerial. Best for showing the reservoirs as water bodies, but
#                  it carries no labels, so the annotation has to do more work.
#   WorldTerrain / WorldShadedRelief  quieter, good if the figure is going in a
#                  report where the river lines must dominate.
BASEMAP = "WorldTopoMap"
# A second, label-only layer drawn over aerial imagery. Ignored unless BASEMAP
# is WorldImagery, where it puts place names back on top. Given as an explicit
# URL because xyzservices does not carry Esri's reference layers -- only the
# nine full basemaps above -- so there is no cx.providers.Esri entry to name.
# Esri tile URLs are {z}/{y}/{x}, not the usual {z}/{x}/{y}.
IMAGERY_LABEL_OVERLAY = True
ESRI_REFERENCE_URL = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
                      "Reference/World_Boundaries_and_Places/MapServer/"
                      "tile/{z}/{y}/{x}")
# Higher = sharper basemap and more tiles fetched. None lets contextily pick
# from the extent. 10-11 suits a basin this size.
ZOOM = 10
BASEMAP_ALPHA = 1.0

# --- what the map covers -----------------------------------------------------
# Padding around the basin, as a fraction of its width/height.
EXTENT_PAD = 0.06
# Clip the drawn flowlines to the basin polygon. NLDI navigation can return a
# little geometry outside the delineated basin at the edges.
CLIP_TO_BASIN = True

# --- NLDI / NWIS -------------------------------------------------------------
# NLDI moved hosts; try them in order and use the first that answers.
NLDI_BASES = [
    "https://api.water.usgs.gov/nldi/linked-data",
    "https://labs.waterdata.usgs.gov/api/nldi/linked-data",
]
NWIS_SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
HTTP_TIMEOUT = 90
# km upstream to navigate. The Cowlitz above Castle Rock is ~ 105 river miles;
# 400 km is comfortably past the headwaters and still bounded.
NAV_DISTANCE_KM = 400

# The outlet the basin and the mainstem are navigated from.
OUTLET_SITE = "14243000"

# --- gages -------------------------------------------------------------------
# (NWIS number, label, sublabel, label offset in points)
#
# The offsets matter more than they look. The Mayfield outflow gage is about
# 800 m below Mayfield Dam, so at basin scale the gage marker and the dam
# marker are the same dot: their labels have to be thrown in opposite
# directions or they overprint. Anything offset far enough to clear gets a
# leader line back to its marker (see LEADER_MIN_POINTS), which is what keeps
# a displaced label from looking like it belongs to the wrong feature.
#
# Keep the sublabels short. The story -- which dam holds the flood, why the
# Toutle matters -- belongs in the figure caption, not printed over the map.
GAGES = [
    ("14243000", "Castle Rock", "USGS 14243000", (20, -6)),
    ("14238000", "Mayfield outflow gage", "USGS 14238000", (-30, 34)),
]
# Seeds the Toutle line. Toutle R. at Tower Rd near Silver Lake, the long-record
# station near the mouth. Set to None to leave the Toutle in the context network.
TOUTLE_SITE = "14242580"

# --- dams and reservoirs -----------------------------------------------------
# NOT NWIS sites, so these are not fetched -- they are typed in. They are only
# used to place a marker and a label; no number in the study depends on them.
# Cross-check against the basemap once and nudge if a marker sits off the dam.
DAMS = [
    ("Mossyrock Dam", 46.5347, -122.4331, "1968", (18, 24)),
    ("Mayfield Dam", 46.5033, -122.5883, "1963", (-18, -44)),
]
# Labels only; the basemap draws the water. (lat, lon, label, rotation)
LAKES = [
    (46.5300, -122.3300, "Riffe Lake", -10),
    (46.5140, -122.5250, "Mayfield Lake", -5),
]
SHOW_DAMS = True
SHOW_LAKE_LABELS = True

# Used ONLY if the NWIS fetch fails, so the figure still has its gages. These
# are approximate; anything drawn from them is stamped on the figure.
FALLBACK_SITES = {
    "14243000": (46.2751, -122.9051),
    "14238000": (46.5028, -122.5992),
    "14242580": (46.3369, -122.8747),
}

# --- styling -----------------------------------------------------------------
C_MAIN = "#1a4f8a"        # Cowlitz mainstem
C_TRIB = "#5b9bd5"        # Toutle
C_NET = "#8fbcdb"         # the rest of the network
C_BASIN = "#22313f"       # basin outline
C_GAGE = "#c0392b"        # gage markers
C_DAM = "#2c3e50"         # dam markers
LW_MAIN, LW_TRIB, LW_NET = 3.2, 2.2, 0.7
BASIN_FILL_ALPHA = 0.07
LABEL_HALO = 2.6          # white outline on text, so labels read over imagery
# A label offset at least this far (in points) gets a leader line back to its
# marker. Below it the label sits close enough to read as attached already.
LEADER_MIN_POINTS = 22.0
TITLE = "Cowlitz River Basin - Castle Rock flow frequency study"
SUBTITLE = ("Flowlines: USGS NHDPlus via NLDI.  Gage locations: USGS NWIS.  "
            "Basin: NLDI delineation above 14243000.")
SHOW_SCALEBAR = True
SCALEBAR_KM = 20
# Scale bar anchor as a fraction of the axes, from the lower left. Nudge it if
# it lands under the Castle Rock label -- both want the bottom-left corner.
SCALEBAR_ANCHOR = (0.055, 0.085)
SHOW_NORTH_ARROW = True

# ----------------------------------------------------------------------------

WGS84 = "EPSG:4326"
WEBM = "EPSG:3857"


def cache_path(name):
    return os.path.join(CACHE_DIR, name)


def ensure_cache_dir():
    if not os.path.isdir(CACHE_DIR):
        os.makedirs(CACHE_DIR)


def http_get(url, params=None):
    """One GET, returning text, or None with the reason printed.

    requests is imported here rather than at module scope so the script can be
    read and its settings edited in an environment that does not have it.
    """
    import requests
    try:
        response = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
        response.raise_for_status()
        return response.text
    except Exception as exc:
        print("      fetch failed: %s" % exc)
        return None


def nldi_get(path, params=None):
    """GET an NLDI path, trying each host in NLDI_BASES."""
    for base in NLDI_BASES:
        url = "%s/%s" % (base, path.lstrip("/"))
        print("   NLDI %s" % url)
        text = http_get(url, params)
        if text:
            try:
                data = json.loads(text)
            except ValueError:
                print("      not JSON, trying the next host")
                continue
            if data and data.get("features"):
                return data
            print("      no features, trying the next host")
    return None


def read_cached_geojson(name):
    path = cache_path(name)
    if not os.path.isfile(path):
        return None
    try:
        frame = gpd.read_file(path)
        return frame if len(frame) else None
    except Exception as exc:
        print("   cache %s unreadable (%s), refetching" % (name, exc))
        return None


def fetch_flowlines(name, site, mode):
    """Navigated flowlines as a GeoDataFrame, from the cache or from NLDI.

    mode is the NLDI navigation code: UM = upstream mainstem, UT = upstream
    with tributaries.
    """
    cached = read_cached_geojson(name)
    if cached is not None:
        print("   %-22s cached (%d features)" % (name, len(cached)))
        return cached
    data = nldi_get("nwissite/USGS-%s/navigation/%s/flowlines" % (site, mode),
                    {"f": "json", "distance": NAV_DISTANCE_KM})
    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None
    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")
    print("   %-22s fetched (%d features)" % (name, len(frame)))
    return frame


def fetch_basin(name, site):
    cached = read_cached_geojson(name)
    if cached is not None:
        print("   %-22s cached" % name)
        return cached
    data = nldi_get("nwissite/USGS-%s/basin" % site, {"f": "json"})
    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None
    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")
    print("   %-22s fetched" % name)
    return frame


def fetch_sites(numbers):
    """Gage coordinates from NWIS, cached. Falls back to FALLBACK_SITES.

    Returns (frame, used_fallback). The flag is carried to the figure rather
    than being swallowed, because a map drawn from typed-in coordinates should
    say so.
    """
    path = cache_path("sites.csv")
    if os.path.isfile(path):
        table = pd.read_csv(path, dtype={"site_no": str})
        has_fallback = bool(table.get("is_fallback", pd.Series([False])).any())
        if set(numbers) <= set(table["site_no"]) and not has_fallback:
            print("   %-22s cached" % "gage coordinates")
            return table, False
        if has_fallback:
            # Never let a fallback coordinate survive in the cache. It was
            # written by a run that could not reach NWIS, and without this the
            # cache check would keep serving typed-in coordinates as though
            # they were published ones long after the network came back.
            print("   %-22s cached copy holds fallbacks -- refetching"
                  % "gage coordinates")

    print("   NWIS %s" % NWIS_SITE_URL)
    text = http_get(NWIS_SITE_URL, {"format": "rdb", "sites": ",".join(numbers),
                                    "siteOutput": "expanded"})
    rows = []
    if text:
        lines = [ln for ln in text.splitlines()
                 if ln and not ln.startswith("#")]
        if len(lines) >= 3:
            header = lines[0].split("\t")
            for line in lines[2:]:                    # line 1 is the format row
                parts = line.split("\t")
                if len(parts) != len(header):
                    continue
                rec = dict(zip(header, parts))
                try:
                    rows.append({"site_no": rec["site_no"].strip(),
                                 "station_nm": rec.get("station_nm", "").strip(),
                                 "lat": float(rec["dec_lat_va"]),
                                 "lon": float(rec["dec_long_va"]),
                                 "is_fallback": False})
                except (KeyError, ValueError):
                    continue

    got = {r["site_no"] for r in rows}
    missing = [n for n in numbers if n not in got]
    for number in missing:
        if number in FALLBACK_SITES:
            lat, lon = FALLBACK_SITES[number]
            rows.append({"site_no": number, "station_nm": "(fallback)",
                         "lat": lat, "lon": lon, "is_fallback": True})
    if missing:
        print("   *** NWIS did not return %s -- using FALLBACK_SITES, which are"
              % ", ".join(missing))
        print("       approximate. The figure is stamped accordingly.")

    table = pd.DataFrame(rows)
    if len(table):
        ensure_cache_dir()
        table.to_csv(path, index=False)
    print("   %-22s %d site(s)" % ("gage coordinates", len(table)))
    return table, bool(len(table)) and bool(table["is_fallback"].any())


def to_web(frame):
    return None if frame is None or not len(frame) else frame.to_crs(WEBM)


def clip_to(frame, basin):
    """Clip to the basin, tolerating the odd invalid ring NLDI can return."""
    if frame is None or basin is None or not CLIP_TO_BASIN:
        return frame
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return gpd.clip(frame, basin)
    except Exception as exc:
        print("   clip skipped (%s)" % exc)
        return frame


def points_frame(table):
    return gpd.GeoDataFrame(
        table.copy(),
        geometry=[Point(x, y) for x, y in zip(table["lon"], table["lat"])],
        crs=WGS84).to_crs(WEBM)


def annotate_leader(ax, xy, text, offset, color, size=9):
    """Label a marker, with a leader line back to it when it sits far off.

    A label nudged a few points off its marker reads as attached. One thrown 30
    points away to dodge a neighbour does not, and on a cluster like Mayfield
    Dam plus its outflow gage an unattached label is worse than no label -- it
    reads as belonging to whichever marker it happens to land nearest.
    """
    far = (offset[0] ** 2 + offset[1] ** 2) ** 0.5 >= LEADER_MIN_POINTS
    ax.annotate(text, xy=xy, xytext=offset, textcoords="offset points",
                ha="left" if offset[0] >= 0 else "right",
                va="bottom" if offset[1] >= 0 else "top",
                zorder=11,
                arrowprops=dict(arrowstyle="-", color=color, lw=0.9,
                                shrinkA=0, shrinkB=6, alpha=0.85)
                if far else None,
                **halo(size, "bold", color))


def halo(size, weight="normal", color="black"):
    """Text kwargs with a white outline, so labels survive a busy basemap."""
    import matplotlib.patheffects as pe
    return dict(fontsize=size, fontweight=weight, color=color,
                path_effects=[pe.withStroke(linewidth=LABEL_HALO,
                                            foreground="white")])


def add_scalebar(ax, length_km=None):
    """Plain scale bar. Web Mercator distorts scale with latitude, so the bar
    length is corrected for the latitude at the middle of the map -- otherwise
    it reads about 45% long at 46 N."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    length_km = SCALEBAR_KM if length_km is None else length_km
    lat = np.degrees(2 * np.arctan(np.exp((y0 + y1) / 2 / 6378137.0))
                     - np.pi / 2)
    metres = length_km * 1000.0 / np.cos(np.radians(lat))
    x_start = x0 + SCALEBAR_ANCHOR[0] * (x1 - x0)
    y_pos = y0 + SCALEBAR_ANCHOR[1] * (y1 - y0)
    ax.plot([x_start, x_start + metres], [y_pos, y_pos], color="black", lw=3,
            solid_capstyle="butt", zorder=12)
    for x in (x_start, x_start + metres):
        ax.plot([x, x], [y_pos, y_pos + 0.012 * (y1 - y0)], color="black",
                lw=3, zorder=12)
    ax.text(x_start + metres / 2, y_pos + 0.018 * (y1 - y0),
            "%d km" % length_km, ha="center", va="bottom", **halo(9, "bold"))


def add_north_arrow(ax):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x = x0 + 0.045 * (x1 - x0)
    y = y1 - 0.135 * (y1 - y0)
    ax.annotate("N", xy=(x, y), xytext=(x, y - 0.055 * (y1 - y0)),
                ha="center", va="center", zorder=12,
                arrowprops=dict(arrowstyle="-|>", color="black", lw=2),
                **halo(12, "bold"))


def main():
    ensure_cache_dir()
    print("=" * 78)
    print("Basemap   : Esri.%s (zoom %s)" % (BASEMAP, ZOOM))
    print("Cache     : %s" % os.path.abspath(CACHE_DIR))
    print("=" * 78)

    wanted = [n for n, _, _, _ in GAGES]
    if TOUTLE_SITE and TOUTLE_SITE not in wanted:
        wanted.append(TOUTLE_SITE)
    sites, used_fallback = fetch_sites(wanted)

    basin = fetch_basin("basin.geojson", OUTLET_SITE)
    network = fetch_flowlines("network_ut.geojson", OUTLET_SITE, "UT")
    mainstem = fetch_flowlines("cowlitz_um.geojson", OUTLET_SITE, "UM")
    toutle = (fetch_flowlines("toutle_um.geojson", TOUTLE_SITE, "UM")
              if TOUTLE_SITE else None)

    if all(x is None for x in (basin, network, mainstem, toutle)):
        raise SystemExit(
            "No geometry was fetched and the cache is empty, so there is "
            "nothing to draw.\nCheck the network, then re-run. The gage "
            "coordinates alone are not a map.")

    basin_w = to_web(basin)
    network_w = clip_to(to_web(network), basin_w)
    mainstem_w = clip_to(to_web(mainstem), basin_w)
    toutle_w = clip_to(to_web(toutle), basin_w)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # --- extent, set before the basemap so contextily fetches the right tiles
    frame_for_extent = basin_w if basin_w is not None else network_w
    minx, miny, maxx, maxy = frame_for_extent.total_bounds
    dx, dy = (maxx - minx) * EXTENT_PAD, (maxy - miny) * EXTENT_PAD
    ax.set_xlim(minx - dx, maxx + dx)
    ax.set_ylim(miny - dy, maxy + dy)

    if basin_w is not None:
        basin_w.plot(ax=ax, facecolor=C_BASIN, alpha=BASIN_FILL_ALPHA,
                     edgecolor="none", zorder=2)
        basin_w.boundary.plot(ax=ax, color=C_BASIN, lw=1.6, ls="--", zorder=3)
    if network_w is not None:
        network_w.plot(ax=ax, color=C_NET, lw=LW_NET, zorder=4)
    if mainstem_w is not None:
        mainstem_w.plot(ax=ax, color=C_MAIN, lw=LW_MAIN, zorder=6)
    if toutle_w is not None:
        toutle_w.plot(ax=ax, color=C_TRIB, lw=LW_TRIB, zorder=5)

    # --- gages
    gage_pts = points_frame(sites[sites["site_no"].isin(
        [n for n, _, _, _ in GAGES])]) if len(sites) else None
    if gage_pts is not None:
        by_site = {r["site_no"]: r for _, r in gage_pts.iterrows()}
        for number, label, sub, offset in GAGES:
            if number not in by_site:
                print("   gage %s has no coordinate -- not drawn" % number)
                continue
            row = by_site[number]
            ax.plot([row.geometry.x], [row.geometry.y], marker="o", ms=11,
                    mfc=C_GAGE, mec="white", mew=1.8, zorder=10)
            annotate_leader(ax, (row.geometry.x, row.geometry.y),
                            "%s\n%s" % (label, sub), offset, C_GAGE)

    # --- dams
    if SHOW_DAMS:
        dam_pts = points_frame(pd.DataFrame(
            [{"lat": la, "lon": lo} for _, la, lo, _, _ in DAMS]))
        for (name, _, _, note, offset), (_, row) in zip(DAMS,
                                                        dam_pts.iterrows()):
            ax.plot([row.geometry.x], [row.geometry.y], marker="s", ms=10,
                    mfc=C_DAM, mec="white", mew=1.6, zorder=10)
            annotate_leader(ax, (row.geometry.x, row.geometry.y),
                            "%s\n%s" % (name, note), offset, C_DAM)

    # --- lake labels
    if SHOW_LAKE_LABELS:
        lake_pts = points_frame(pd.DataFrame(
            [{"lat": la, "lon": lo} for la, lo, _, _ in LAKES]))
        for (_, _, label, rot), (_, row) in zip(LAKES, lake_pts.iterrows()):
            ax.text(row.geometry.x, row.geometry.y, label, rotation=rot,
                    ha="center", va="center", style="italic", zorder=11,
                    **halo(9.5, "normal", "#12507a"))

    # --- basemap last, so its extent matches the data already drawn
    try:
        if not hasattr(cx.providers.Esri, BASEMAP):
            raise SystemExit(
                "BASEMAP is '%s', which this xyzservices does not have.\n"
                "Available Esri basemaps: %s"
                % (BASEMAP, ", ".join(sorted(cx.providers.Esri.keys()))))
        source = getattr(cx.providers.Esri, BASEMAP)
        cx.add_basemap(ax, source=source, zoom=ZOOM, alpha=BASEMAP_ALPHA,
                       attribution_size=6)
        if BASEMAP == "WorldImagery" and IMAGERY_LABEL_OVERLAY:
            cx.add_basemap(ax, source=ESRI_REFERENCE_URL, zoom=ZOOM,
                           attribution=False)
    except SystemExit:
        raise
    except Exception as exc:
        print("\n   *** BASEMAP NOT DRAWN: %s" % exc)
        print("   The vectors are still correct; only the backdrop is missing.")
        print("   Usually the tile host is unreachable. contextily caches "
              "tiles, so a")
        print("   later run on a working connection will fill it in.")
        ax.set_facecolor("#eef2f5")

    ax.set_xlim(minx - dx, maxx + dx)
    ax.set_ylim(miny - dy, maxy + dy)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#33414d")
        spine.set_linewidth(1.0)

    if SHOW_SCALEBAR:
        add_scalebar(ax)
    if SHOW_NORTH_ARROW:
        add_north_arrow(ax)

    handles = [
        Line2D([], [], color=C_MAIN, lw=LW_MAIN, label="Cowlitz River (mainstem)"),
        Line2D([], [], color=C_TRIB, lw=LW_TRIB,
               label="Toutle River (unregulated tributary)"),
        Line2D([], [], color=C_NET, lw=1.4, label="Other tributaries"),
        Line2D([], [], color=C_GAGE, marker="o", ms=9, mfc=C_GAGE, mec="white",
               ls="none", label="USGS streamgage"),
        Line2D([], [], color=C_DAM, marker="s", ms=9, mfc=C_DAM, mec="white",
               ls="none", label="Dam"),
        Patch(facecolor=C_BASIN, alpha=0.18, edgecolor=C_BASIN, ls="--",
              label="Basin above Castle Rock"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=9, framealpha=0.93)

    ax.set_title(TITLE, fontsize=13.5, fontweight="bold", pad=10)
    note = SUBTITLE
    if used_fallback:
        note += ("\nWARNING: one or more gage positions came from "
                 "FALLBACK_SITES (approximate), not from NWIS.")
    ax.text(0.5, -0.035, note, transform=ax.transAxes, ha="center", va="top",
            fontsize=7.5, color="#c0392b" if used_fallback else "#4a5568")

    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("\nWrote %s  (%.1f x %.1f in at %d dpi)"
          % (os.path.abspath(OUT_PNG), FIGSIZE[0], FIGSIZE[1], DPI))
    if used_fallback:
        print("*** At least one gage used FALLBACK_SITES. Delete %s and re-run"
              % cache_path("sites.csv"))
        print("    with a working connection to get the published coordinates.")


main()
