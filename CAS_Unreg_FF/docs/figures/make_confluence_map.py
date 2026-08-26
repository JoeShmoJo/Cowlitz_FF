#make_confluence_map.py
# -*- coding: utf-8 -*-
"""
Cowlitz River location map extended DOWNSTREAM to the Coweeman confluence,
showing Arkansas Creek, Ostrander Creek and the Coweeman River.

A COPY of make_basin_map.py. The original still makes the Castle Rock basin
map for the unregulated memo and is untouched; this one supports Section 8 of
the combined memo, where the regulated curve is carried from the gage down to
the Coweeman confluence. Separate output, separate cache -- running one does
not disturb the other.

WHAT IS DIFFERENT FROM THE ORIGINAL
    - the Cowlitz is navigated DOWNSTREAM from Castle Rock as well as up, so
      the mainstem reaches past all three local confluences
    - Arkansas Creek, Ostrander Creek and the Coweeman River are drawn and
      labelled exactly like the Toutle and the Tilton: one colour, one
      weight, names written along their own lines. On this map they are all
      the same kind of feature.
    - ONE watershed boundary, delineated just BELOW the Coweeman confluence
      so it contains every creek drawn. The Castle Rock basin alone would
      leave all three creeks outside the boundary, which reads as an error.
    - CLIP_TO_BASIN is OFF. The clip is against the basin ABOVE the gage in
      the original; here the added features are below it, and clipping would
      silently delete the entire point of this figure.

SEEDING THE UNGAGED CREEKS
    The original seeds every tributary from an NWIS gage. Arkansas and
    Ostrander have no gage, which is the whole reason Section 8 estimates
    their contribution by drainage area. So this version can also seed from a
    COMID, or from a coordinate that NLDI snaps to the nearest flowline:

        {"site":  "14245000"}      an NWIS station
        {"comid": "23735691"}      an NHDPlus reach
        {"point": (lat, lon)}      snapped to the nearest reach

    THE TWO CREEK COORDINATES BELOW ARE APPROXIMATE SEEDS, NOT SURVEYED
    POINTS. They only have to land near the right creek; NLDI snaps to the
    nearest flowline from there. The script PRINTS the comid it resolved and
    the distance it moved, so an obviously wrong snap is visible in the log.
    Check the drawn lines against the basemap on the first run and nudge the
    coordinates if a creek is wrong. Nothing in any result depends on them --
    this is a figure.

LABELS ARE TIGHT DOWN HERE
    The three confluences are within about 20 km of river and the creeks are
    short, so their names crowd each other and the mainstem. label_frac slides
    a name along its own line and is the dial to turn; which END it counts
    from is not predictable, because linemerge does not preserve direction, so
    treat it as a dial rather than as "0 is the mouth". Expect to adjust the
    three of them once against the real basemap.

RUNNING IT
    Same environment as make_basin_map.py -- environment.yml beside this file.

        conda env create -f environment.yml
        conda activate cowlitz-map
        python make_confluence_map.py
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
from shapely.geometry import Point

# ----------------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------------
OUT_PNG = r"confluence_map.png"
# Separate cache: this map fetches different geometry and must not
# overwrite what make_basin_map.py has already stored.
CACHE_DIR = r"mapdata_confluence"
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
# MUST stay False here. Everything this map adds -- the downstream mainstem
# and all three creeks -- is BELOW the Castle Rock gage and so outside the
# basin delineated above it. Clipping to that basin erases it all.
CLIP_TO_BASIN = False

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

# The outlet the basin and the UPSTREAM mainstem are navigated from.
OUTLET_SITE = "14243000"

# How far to carry the mainstem BELOW the gage. Castle Rock to the Coweeman
# confluence is roughly 20 km of river; 45 km clears it with margin without
# running all the way to the Columbia.
DOWNSTREAM_KM = 45

# Local tributaries below the gage, drawn and labelled. seed is one of:
#     {"site": "14245000"}   an NWIS station
#     {"comid": "23735691"}  an NHDPlus reach
#     {"point": (lat, lon)}  snapped to the nearest reach by NLDI
#
# THE TWO POINT SEEDS ARE APPROXIMATE. Arkansas and Ostrander have no gage --
# that is exactly why Section 8 estimates them by drainage area -- so there is
# no station to seed from. The coordinates only have to land near the right
# creek; the run prints the comid resolved and how far the snap moved, and the
# drawn line should be checked against the basemap the first time. If a creek
# comes out wrong, nudge its coordinate. No result depends on this.
LOCAL_TRIBUTARIES = [
    {"label": "Coweeman River", "seed": {"site": "14245000"},
     "up_km": 90, "down_km": 30, "label_frac": 0.45},
    {"label": "Ostrander Creek", "seed": {"point": (46.2085, -122.8790)},
     "up_km": 40, "down_km": 12, "label_frac": 0.40},
    {"label": "Arkansas Creek", "seed": {"point": (46.2430, -122.8880)},
     "up_km": 40, "down_km": 12, "label_frac": 0.40},
]
# The creeks are drawn and named exactly like the Toutle and the Tilton --
# same colour, same weight, names written along their own lines. On this map
# they are all the same kind of feature and styling them apart would imply a
# distinction that is not there. Crowding is handled with label_frac, which
# slides a name along its own line.

# ONE watershed boundary: the Cowlitz basin, delineated at a point just BELOW
# the Coweeman confluence so it contains every creek this map draws. The
# Castle Rock basin alone would leave all three creeks outside the boundary,
# which reads as an error even though it is technically correct.
#
# Approximate seed, snapped to the nearest flowline like the creeks. The run
# prints the resulting area -- expect roughly 2,480 sq mi against the 2,476
# Section 8 uses at the confluence. A wildly different number means the seed
# landed on the wrong stream; nudge it and delete the cached basin.
BASIN_SEED = {"point": (46.1300, -122.9080)}
BASIN_EXPECT_SQ_MI = 2476.0

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
    ("14243000", "Castle Rock", "USGS 14243000", (22, -10)),
    ("14238000", "Mayfield outflow gage", "USGS 14238000", (-30, 34)),
    # The Coweeman gage. Its record is what Section 8's drainage-area ratio
    # was tested against, so it earns a marker on this map.
    ("14245000", "Coweeman R nr Kelso", "USGS 14245000", (26, -20)),
]
# --- named tributaries -------------------------------------------------------
# Drawn in one style, heavier than the context network, and labelled along the
# line. These are the tributaries worth picking out: the Toutle joins BELOW
# both dams, so it reaches Castle Rock unregulated, and the Tilton joins
# Mayfield Lake, so it is caught by the projects.
#
# Each is seeded from an NWIS gage and drawn by navigating BOTH ways from it:
#   UM  upstream to the headwaters
#   DM  downstream to the confluence
# Upstream alone is what left the Toutle hanging in mid-air, stopping at its
# gage instead of reaching the Cowlitz.
#
# DM does not stop at the confluence -- it carries on down the Cowlitz mainstem
# to Castle Rock. Those reaches are removed by dropping any comid that is also
# in the Cowlitz mainstem, which is exact (comid identity, not geometry
# matching) and leaves the tributary ending precisely where it meets the river.
#
# site: an NWIS number, or None to look one up in the basin by name_match. A
# lookup prints every candidate it found, so a river with gages on several
# forks can be pinned to the right one by filling in site.
#
# label_frac: position along the merged line, 0 to 1. Which END that starts
# from is NOT predictable -- shapely's linemerge does not preserve a direction
# -- so treat it as a dial to turn until the label sits somewhere sensible,
# not as "0 is the mouth". The Toutle is set away from the middle on purpose:
# its confluence is close to Castle Rock, so a label near the mouth lands on
# top of the Castle Rock gage label.
TRIBUTARIES = [
    {"label": "Toutle River", "site": "14242580", "name_match": "TOUTLE",
     "label_frac": 0.62},
    {"label": "Tilton River", "site": None, "name_match": "TILTON",
     "label_frac": 0.45},
]
# How far downstream to navigate from a tributary gage. Only has to reach the
# confluence; the mainstem reaches beyond it are dropped anyway.
TRIB_DOWNSTREAM_KM = 60
# Label the Cowlitz mainstem along the line as well. With the tributaries named
# on the map there is no legend left to say which line is which.
MAINSTEM_LABEL = "Cowlitz River"
MAINSTEM_LABEL_FRAC = 0.62
# Rotate river labels to follow the line. False leaves them horizontal.
ROTATE_RIVER_LABELS = True

# --- dams and reservoirs -----------------------------------------------------
# NOT NWIS sites, so these are not fetched -- they are typed in. They are only
# used to place a marker and a label; no number in the study depends on them.
# Cross-check against the basemap once and nudge if a marker sits off the dam.
DAMS = [
    ("Mossyrock Dam", 46.5347, -122.4331, "1968", (0, -40)),
    ("Mayfield Dam", 46.5033, -122.5883, "1963", (-10, -30)),
]
# Labels only; the basemap draws the water. (lat, lon, label, rotation)
LAKES = [
    (46.5300, -122.3300, "Riffe Lake", -10),
    (46.5140, -122.5250, "Mayfield Lake", -5),
]
SHOW_DAMS = True
SHOW_LAKE_LABELS = False

# Used ONLY if the NWIS fetch fails, so the figure still has its gages. These
# are approximate; anything drawn from them is stamped on the figure.
# NOTE: no entry for 14245000. A fallback coordinate has to be right to be
# worth having, and an invented one would place the Coweeman gage somewhere
# plausible-looking and wrong. If the NWIS fetch fails, that marker is simply
# not drawn and the run says so.
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
TITLE = ("Lower Cowlitz River - Castle Rock gage to the Coweeman confluence")
SUBTITLE = ("Flowlines: USGS NHDPlus via NLDI.  Gage locations: USGS NWIS.  "
            "Basin: NLDI delineation below the Coweeman confluence.")
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


def fetch_flowlines(name, site, mode, distance=None):
    """Navigated flowlines as a GeoDataFrame, from the cache or from NLDI.

    mode is the NLDI navigation code: UM = upstream mainstem, UT = upstream
    with tributaries, DM = downstream mainstem.
    """
    cached = read_cached_geojson(name)
    if cached is not None:
        print("   %-22s cached (%d features)" % (name, len(cached)))
        return cached
    data = nldi_get("nwissite/USGS-%s/navigation/%s/flowlines" % (site, mode),
                    {"f": "json",
                     "distance": NAV_DISTANCE_KM if distance is None else distance})
    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None
    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")
    print("   %-22s fetched (%d features)" % (name, len(frame)))
    return frame


def comid_set(frame):
    """The NHDPlus comids in a flowline frame, as strings.

    NLDI has used more than one spelling for this property across versions, so
    take whichever is present rather than assuming.
    """
    if frame is None or not len(frame):
        return set()
    for column in ("nhdplus_comid", "comid", "COMID", "nhdplusComid"):
        if column in frame.columns:
            return {str(v) for v in frame[column].dropna()}
    return set()


def drop_comids(frame, unwanted):
    """Remove flowlines whose comid is in `unwanted`."""
    if frame is None or not len(frame) or not unwanted:
        return frame
    for column in ("nhdplus_comid", "comid", "COMID", "nhdplusComid"):
        if column in frame.columns:
            keep = ~frame[column].astype(str).isin(unwanted)
            return frame[keep].copy()
    return frame


def find_site_by_name(pattern, cache_name="upstream_sites.geojson"):
    """An NWIS site in the basin whose name contains `pattern`.

    Used when a tributary has no site number filled in. Every candidate is
    printed, because a river with gages on several forks will match more than
    one and only the operator can say which fork should draw the line.
    """
    frame = read_cached_geojson(cache_name)
    if frame is None:
        data = nldi_get("nwissite/USGS-%s/navigation/UT/nwissite" % OUTLET_SITE,
                        {"f": "json", "distance": NAV_DISTANCE_KM})
        if not data:
            print("   site lookup for '%s' UNAVAILABLE" % pattern)
            return None
        frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
        ensure_cache_dir()
        frame.to_file(cache_path(cache_name), driver="GeoJSON")

    name_col = next((c for c in ("name", "NAME", "station_nm")
                     if c in frame.columns), None)
    id_col = next((c for c in ("identifier", "identifie", "ID")
                   if c in frame.columns), None)
    if not name_col or not id_col:
        print("   site lookup: unexpected columns %s" % list(frame.columns))
        return None

    hits = frame[frame[name_col].astype(str).str.upper().str.contains(
        pattern.upper(), na=False)]
    if not len(hits):
        print("   site lookup: nothing in the basin matches '%s'" % pattern)
        return None
    numbers = []
    for _, row in hits.iterrows():
        number = str(row[id_col]).replace("USGS-", "").strip()
        numbers.append(number)
        print("      candidate %-12s %s" % (number, row[name_col]))
    if len(numbers) > 1:
        print("      using %s -- set 'site' in TRIBUTARIES to pin another"
              % numbers[0])
    return numbers[0]


def fetch_tributary(trib, mainstem_comids):
    """One named tributary, headwaters to confluence.

    UM gives the upstream half. DM gives the downstream half but runs on past
    the confluence and down the Cowlitz, so the mainstem comids are subtracted
    -- exact, and it leaves the line ending where the rivers actually meet.
    """
    label = trib["label"]
    stem = label.split()[0].lower()
    up_name, down_name = "%s_um.geojson" % stem, "%s_dm.geojson" % stem

    # Cache first, seed second. The seed is only needed to FETCH; once the
    # reaches are on disk the river can be drawn with no network and no gage
    # lookup, which is the whole point of the cache.
    up, down = read_cached_geojson(up_name), read_cached_geojson(down_name)
    site = trib.get("site")
    if up is None or down is None:
        if not site:
            print("   %s: no site number, looking one up" % label)
            site = find_site_by_name(trib["name_match"])
        if not site:
            print("   %s NOT DRAWN (no seed gage and nothing cached)" % label)
            return None
        up = fetch_flowlines(up_name, site, "UM")
        down = fetch_flowlines(down_name, site, "DM",
                               distance=TRIB_DOWNSTREAM_KM)
    down = drop_comids(down, mainstem_comids)
    parts = [p for p in (up, down) if p is not None and len(p)]
    if not parts:
        print("   %s NOT DRAWN (no flowlines)" % label)
        return None
    frame = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=WGS84)
    print("   %-22s %d features (seed %s)"
          % (label, len(frame), site or "cached"))
    return frame


def comid_at_point(lat, lon):
    """The NHDPlus comid nearest a coordinate, via NLDI's position lookup.

    NLDI snaps to the nearest flowline, so an approximate coordinate is
    enough. What is printed is the point that came BACK, which is the check
    that matters: if the snap moved a long way, the seed was in the wrong
    place and the creek drawn from it will be the wrong creek.
    """
    data = nldi_get("comid/position", {"coords": "POINT(%f %f)" % (lon, lat)})
    if not data or not data.get("features"):
        print("      position lookup failed for (%.4f, %.4f)" % (lat, lon))
        return None
    feature = data["features"][0]
    comid = None
    for key in ("comid", "nhdplus_comid", "identifier"):
        if feature.get("properties", {}).get(key):
            comid = str(feature["properties"][key])
            break
    if comid is None:
        print("      position lookup returned no comid")
        return None
    try:
        coords = np.array(feature["geometry"]["coordinates"], dtype=float)
        coords = coords.reshape(-1, coords.shape[-1])[:, :2]
        moved = np.min(np.hypot((coords[:, 0] - lon) * np.cos(np.radians(lat)),
                                coords[:, 1] - lat)) * 111.0
        print("      seed (%.4f, %.4f) -> comid %s, snapped %.2f km"
              % (lat, lon, comid, moved))
    except Exception:                                    # noqa: BLE001
        print("      seed (%.4f, %.4f) -> comid %s" % (lat, lon, comid))
    return comid


def seed_path(seed):
    """NLDI path prefix for a seed dict, resolving a point to a comid."""
    if seed.get("site"):
        return "nwissite/USGS-%s" % seed["site"], "site %s" % seed["site"]
    if seed.get("comid"):
        return "comid/%s" % seed["comid"], "comid %s" % seed["comid"]
    if seed.get("point"):
        comid = comid_at_point(*seed["point"])
        if comid is None:
            return None, None
        return "comid/%s" % comid, "comid %s (from point)" % comid
    return None, None


def fetch_flowlines_seeded(name, prefix, mode, distance):
    """Flowlines from any NLDI seed prefix, cached like fetch_flowlines."""
    cached = read_cached_geojson(name)
    if cached is not None:
        print("   %-22s cached (%d features)" % (name, len(cached)))
        return cached
    if prefix is None:
        return None
    data = nldi_get("%s/navigation/%s/flowlines" % (prefix, mode),
                    {"f": "json", "distance": distance})
    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None
    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")
    print("   %-22s fetched (%d features)" % (name, len(frame)))
    return frame


def fetch_local_tributary(trib, mainstem_comids):
    """One creek below the gage, headwaters to confluence.

    Same UM + DM idea as fetch_tributary, but seeded from a site, a comid or a
    snapped coordinate, and with per-creek navigation distances because these
    are small streams next to a 105-mile mainstem.
    """
    label = trib["label"]
    stem = label.split()[0].lower()
    up_name, down_name = "local_%s_um.geojson" % stem, "local_%s_dm.geojson" % stem

    up, down = read_cached_geojson(up_name), read_cached_geojson(down_name)
    note = "cached"
    if up is None or down is None:
        prefix, note = seed_path(trib["seed"])
        if prefix is None:
            print("   %s NOT DRAWN (seed did not resolve)" % label)
            return None
        up = fetch_flowlines_seeded(up_name, prefix, "UM", trib.get("up_km", 40))
        down = fetch_flowlines_seeded(down_name, prefix, "DM",
                                      trib.get("down_km", 15))
    # DM runs on down the Cowlitz past the confluence; drop the shared reaches
    # so the creek ends exactly where it meets the river.
    down = drop_comids(down, mainstem_comids)
    parts = [f for f in (up, down) if f is not None and len(f)]
    if not parts:
        print("   %s NOT DRAWN (no flowlines)" % label)
        return None
    frame = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=WGS84)
    print("   %-22s %d features (seed %s)" % (label, len(frame), note))
    return frame


def fetch_basin_seeded(name, seed, fallback_site):
    """The one watershed boundary, delineated from a seed below the confluence.

    Falls back to the gage's own basin if the seed will not resolve, so a bad
    coordinate degrades to the Castle Rock boundary rather than to no boundary
    at all. Either way the area is printed, which is the check that matters.
    """
    cached = read_cached_geojson(name)
    if cached is None:
        prefix, note = seed_path(seed)
        data = nldi_get("%s/basin" % prefix, {"f": "json"}) if prefix else None
        if not data:
            print("   basin from seed UNAVAILABLE -- falling back to the "
                  "basin above gage %s" % fallback_site)
            return fetch_basin(name, fallback_site)
        cached = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
        ensure_cache_dir()
        cached.to_file(cache_path(name), driver="GeoJSON")
        print("   %-22s fetched (seed %s)" % (name, note))
    else:
        print("   %-22s cached" % name)
    try:
        area = cached.to_crs("EPSG:5070").area.sum() / 2.589988e6
        flag = ("" if abs(area - BASIN_EXPECT_SQ_MI) / BASIN_EXPECT_SQ_MI < 0.15
                else "   <-- CHECK: expected about %.0f" % BASIN_EXPECT_SQ_MI)
        print("   %-22s %.0f sq mi%s" % ("basin area", area, flag))
    except Exception:                                     # noqa: BLE001
        pass
    return cached


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


def longest_strand(frame):
    """The river itself, out of NLDI's unordered pile of short reaches.

    Merge them, then take the longest resulting strand -- not whichever reach
    happened to come first. Returns None if there is nothing usable.
    """
    from shapely.ops import linemerge
    if frame is None or not len(frame):
        return None
    merged = linemerge(list(frame.geometry))
    line = (max(merged.geoms, key=lambda g: g.length)
            if merged.geom_type == "MultiLineString" else merged)
    return line if line.length > 0 else None


def label_along_line(ax, frame, text, frac, color, size=10):
    """Write a river name along its own line, angled to follow it.

    The pieces come back from NLDI as many short reaches in no useful order, so
    they are merged and the longest resulting strand is used -- that is the
    river itself rather than whichever reach happened to be first. The angle is
    taken from a chord either side of the label point, long enough not to pick
    up the wiggle of a single reach.
    """
    line = longest_strand(frame)
    if line is None:
        return

    point = line.interpolate(frac, normalized=True)
    rotation = 0.0
    if ROTATE_RIVER_LABELS:
        step = 0.04
        before = line.interpolate(max(frac - step, 0.0), normalized=True)
        after = line.interpolate(min(frac + step, 1.0), normalized=True)
        rotation = np.degrees(np.arctan2(after.y - before.y, after.x - before.x))
        # keep text upright: never let a label read upside down
        if rotation > 90:
            rotation -= 180
        elif rotation < -90:
            rotation += 180
    ax.text(point.x, point.y, text, rotation=rotation, rotation_mode="anchor",
            ha="center", va="bottom", zorder=11, style="italic",
            **halo(size, "bold", color))


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

    sites, used_fallback = fetch_sites([n for n, _, _, _ in GAGES])

    basin = fetch_basin_seeded("basin.geojson", BASIN_SEED, OUTLET_SITE)
    network = fetch_flowlines("network_ut.geojson", OUTLET_SITE, "UT")
    mainstem = fetch_flowlines("cowlitz_um.geojson", OUTLET_SITE, "UM")

    # The mainstem comids are what trims each tributary at its confluence, so
    # they have to be in hand before the tributaries are built.
    mainstem_comids = comid_set(mainstem)
    if mainstem is not None and not mainstem_comids:
        print("   NOTE: no comid column on the mainstem, so the tributaries "
              "cannot be")
        print("         trimmed at their confluences and will run on down the "
              "Cowlitz.")
    tribs = [(t, fetch_tributary(t, mainstem_comids)) for t in TRIBUTARIES]

    # Everything below the gage. The downstream mainstem must be fetched
    # BEFORE the creeks, because its comids are what trim each creek at its
    # confluence -- the upstream mainstem alone does not reach down here.
    downstream = fetch_flowlines_seeded(
        "cowlitz_dm.geojson", "nwissite/USGS-%s" % OUTLET_SITE,
        "DM", DOWNSTREAM_KM)
    local_trim = mainstem_comids | comid_set(downstream)
    locals_ = [(t, fetch_local_tributary(t, local_trim))
               for t in LOCAL_TRIBUTARIES]

    # Include everything this map adds, or a run that fetched the creeks but
    # not the upstream basin would bail out with a map it could have drawn.
    if all(x is None for x in [basin, network, mainstem, downstream]
           + [f for _, f in tribs]
           + [f for _, f in locals_]):
        raise SystemExit(
            "No geometry was fetched and the cache is empty, so there is "
            "nothing to draw.\nCheck the network, then re-run. The gage "
            "coordinates alone are not a map.")

    basin_w = to_web(basin)
    network_w = clip_to(to_web(network), basin_w)
    mainstem_w = clip_to(to_web(mainstem), basin_w)
    tribs_w = [(t, clip_to(to_web(f), basin_w)) for t, f in tribs]
    # NOT clipped -- these lie below the gage, outside the basin above it.
    downstream_w = to_web(downstream)
    locals_w = [(t, to_web(f)) for t, f in locals_]

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # --- extent, set before the basemap so contextily fetches the right tiles
    # The extent has to span the basin above the gage AND everything added
    # below it, or the whole reason for this map falls off the bottom edge.
    extent_parts = [f for f in ([basin_w, network_w, downstream_w]
                                + [f for _, f in locals_w])
                    if f is not None and len(f)]
    if not extent_parts:
        raise SystemExit("nothing to set an extent from")
    bounds = np.array([f.total_bounds for f in extent_parts])
    minx, miny = bounds[:, 0].min(), bounds[:, 1].min()
    maxx, maxy = bounds[:, 2].max(), bounds[:, 3].max()
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
    if downstream_w is not None and len(downstream_w):
        downstream_w.plot(ax=ax, color=C_MAIN, lw=LW_MAIN, zorder=6)
    # same style as the Toutle and the Tilton -- one kind of feature
    for _, frame in locals_w:
        if frame is not None and len(frame):
            frame.plot(ax=ax, color=C_TRIB, lw=LW_TRIB, zorder=5)
    # every named tributary in ONE style -- they are the same kind of thing on
    # this map, and colouring them differently would imply a distinction that
    # is not there
    for _, frame in tribs_w:
        if frame is not None and len(frame):
            frame.plot(ax=ax, color=C_TRIB, lw=LW_TRIB, zorder=5)

    # --- river names, written along the lines instead of into a legend
    for trib, frame in tribs_w:
        label_along_line(ax, frame, trib["label"], trib.get("label_frac", 0.45),
                         C_TRIB)
    for trib, frame in locals_w:
        label_along_line(ax, frame, trib["label"],
                         trib.get("label_frac", 0.45), C_TRIB)
    if MAINSTEM_LABEL and mainstem_w is not None:
        label_along_line(ax, mainstem_w, MAINSTEM_LABEL, MAINSTEM_LABEL_FRAC,
                         C_MAIN, size=11)

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

    # No legend. Every river carries its name along its own line and every gage
    # and dam is labelled at its marker, so a legend would only repeat them --
    # and it covered a corner of the basin to do it.

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
