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
CACHE_DIR = r"mapdata_confluence"
DPI = 300
FIGSIZE = (12.0, 8.0)

# --- basemap -----------------------------------------------------------------
BASEMAP = "WorldTopoMap"
IMAGERY_LABEL_OVERLAY = True
ESRI_REFERENCE_URL = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
                      "Reference/World_Boundaries_and_Places/MapServer/"
                      "tile/{z}/{y}/{x}")
ZOOM = 10
BASEMAP_ALPHA = 1.0

# --- what the map covers -----------------------------------------------------
EXTENT_PAD = 0.06
CLIP_TO_BASIN = False

# --- NLDI / NWIS -------------------------------------------------------------
NLDI_BASES = [
    "https://api.water.usgs.gov/nldi/linked-data",
    "https://labs.waterdata.usgs.gov/api/nldi/linked-data",
]
NWIS_SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
HTTP_TIMEOUT = 90
NAV_DISTANCE_KM = 400

OUTLET_SITE = "14243000"
DOWNSTREAM_KM = 45

LOCAL_TRIBUTARIES = [
    {"label": "Coweeman River", "seed": {"site": "14245000"},
     "up_km": 90, "down_km": 30, "label_frac": 0.45,
     "label_offset": (20, -10),
     "label_rotation": 0},

    {"label": "Ostrander Creek", "seed": {"point": (46.2085, -122.8790)},
     "up_km": 40, "down_km": 12, "label_frac": 0.30,
     "label_offset": (28, -15),
     "label_rotation": 0},

    {"label": "Arkansas Creek", "seed": {"point": (46.2430, -122.8880)},
     "up_km": 40, "down_km": 12, "label_frac": 0.50,
     "label_offset": (45, 0),
     "label_rotation": 0},
]

BASIN_SEED = {"point": (46.09414, -122.93052)}
BASIN_EXPECT_SQ_MI = 2480.0

# --- gages -------------------------------------------------------------------
GAGES = [
    ("14243000", "Castle Rock", "USGS 14243000", (-10, 0)),
    ("14238000", "Mayfield outflow gage", "USGS 14238000", (-30, 34)),
    ("14245000", "Coweeman R nr Kelso", "USGS 14245000", (26, -20)),
]

# --- named tributaries -------------------------------------------------------
TRIBUTARIES = [
    {"label": "Toutle River", "site": "14242580", "name_match": "TOUTLE",
     "label_frac": 0.62},
    {"label": "Tilton River", "site": None, "name_match": "TILTON",
     "label_frac": 0.45},
]

TRIB_DOWNSTREAM_KM = 60

MAINSTEM_LABEL = "Cowlitz River"
MAINSTEM_LABEL_FRAC = 0.62
ROTATE_RIVER_LABELS = True

# --- dams and reservoirs -----------------------------------------------------
DAMS = [
    ("Mossyrock Dam", 46.5347, -122.4331, "", (0, -40)),
    ("Mayfield Dam", 46.5033, -122.5883, "", (-10, -30)),
]

LAKES = [
    (46.5300, -122.3300, "Riffe Lake", -10),
    (46.5140, -122.5250, "Mayfield Lake", -5),
]

SHOW_DAMS = True
SHOW_LAKE_LABELS = False

FALLBACK_SITES = {
    "14243000": (46.2751, -122.9051),
    "14238000": (46.5028, -122.5992),
    "14242580": (46.3369, -122.8747),
}

# --- styling -----------------------------------------------------------------
C_MAIN = "#1a4f8a"
C_TRIB = "#5b9bd5"
C_NET = "#8fbcdb"
C_BASIN = "#22313f"
C_GAGE = "#c0392b"
C_DAM = "#2c3e50"

LW_MAIN, LW_TRIB, LW_NET = 3.2, 2.2, 0.7
BASIN_FILL_ALPHA = 0.07
LABEL_HALO = 2.6
LEADER_MIN_POINTS = 22.0

TITLE = ("Cowlitz River Basin")
SUBTITLE = ("Flowlines: USGS NHDPlus via NLDI.  Gage locations: USGS NWIS.  "
            "Basin: NLDI delineation below the Coweeman confluence.")

SHOW_SCALEBAR = True
SCALEBAR_MI = 10
SCALEBAR_ANCHOR = (0.02, 0.085)
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
    import requests
    try:
        response = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
        response.raise_for_status()
        return response.text
    except Exception as exc:
        print("      fetch failed: %s" % exc)
        return None


def nldi_get(path, params=None):
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
    cached = read_cached_geojson(name)
    if cached is not None:
        print("   %-22s cached (%d features)" % (name, len(cached)))
        return cached

    data = nldi_get(
        "nwissite/USGS-%s/navigation/%s/flowlines" % (site, mode),
        {
            "f": "json",
            "distance": NAV_DISTANCE_KM if distance is None else distance
        }
    )

    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None

    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")
    print("   %-22s fetched (%d features)" % (name, len(frame)))
    return frame


def comid_set(frame):
    if frame is None or not len(frame):
        return set()

    for column in ("nhdplus_comid", "comid", "COMID", "nhdplusComid"):
        if column in frame.columns:
            return {str(v) for v in frame[column].dropna()}

    return set()


def drop_comids(frame, unwanted):
    if frame is None or not len(frame) or not unwanted:
        return frame

    for column in ("nhdplus_comid", "comid", "COMID", "nhdplusComid"):
        if column in frame.columns:
            keep = ~frame[column].astype(str).isin(unwanted)
            return frame[keep].copy()

    return frame


def find_site_by_name(pattern, cache_name="upstream_sites.geojson"):
    frame = read_cached_geojson(cache_name)

    if frame is None:
        data = nldi_get(
            "nwissite/USGS-%s/navigation/UT/nwissite" % OUTLET_SITE,
            {"f": "json", "distance": NAV_DISTANCE_KM}
        )

        if not data:
            print("   site lookup for '%s' UNAVAILABLE" % pattern)
            return None

        frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)
        ensure_cache_dir()
        frame.to_file(cache_path(cache_name), driver="GeoJSON")

    name_col = next(
        (c for c in ("name", "NAME", "station_nm") if c in frame.columns),
        None
    )
    id_col = next(
        (c for c in ("identifier", "identifie", "ID") if c in frame.columns),
        None
    )

    if not name_col or not id_col:
        print("   site lookup: unexpected columns %s" % list(frame.columns))
        return None

    hits = frame[
        frame[name_col].astype(str).str.upper().str.contains(
            pattern.upper(), na=False
        )
    ]

    if not len(hits):
        print("   site lookup: nothing in the basin matches '%s'" % pattern)
        return None

    numbers = []

    for _, row in hits.iterrows():
        number = str(row[id_col]).replace("USGS-", "").strip()
        numbers.append(number)
        print("      candidate %-12s %s" % (number, row[name_col]))

    if len(numbers) > 1:
        print(
            "      using %s -- set 'site' in TRIBUTARIES to pin another"
            % numbers[0]
        )

    return numbers[0]


def fetch_tributary(trib, mainstem_comids):
    label = trib["label"]
    stem = label.split()[0].lower()
    up_name = "%s_um.geojson" % stem
    down_name = "%s_dm.geojson" % stem

    up = read_cached_geojson(up_name)
    down = read_cached_geojson(down_name)

    site = trib.get("site")

    if up is None or down is None:
        if not site:
            print("   %s: no site number, looking one up" % label)
            site = find_site_by_name(trib["name_match"])

        if not site:
            print("   %s NOT DRAWN (no seed gage and nothing cached)" % label)
            return None

        up = fetch_flowlines(up_name, site, "UM")
        down = fetch_flowlines(
            down_name,
            site,
            "DM",
            distance=TRIB_DOWNSTREAM_KM
        )

    down = drop_comids(down, mainstem_comids)

    parts = [
        p for p in (up, down)
        if p is not None and len(p)
    ]

    if not parts:
        print("   %s NOT DRAWN (no flowlines)" % label)
        return None

    frame = gpd.GeoDataFrame(
        pd.concat(parts, ignore_index=True),
        crs=WGS84
    )

    print(
        "   %-22s %d features (seed %s)"
        % (label, len(frame), site or "cached")
    )

    return frame


def comid_at_point(lat, lon):
    data = nldi_get(
        "comid/position",
        {"coords": "POINT(%f %f)" % (lon, lat)}
    )

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
        coords = np.array(
            feature["geometry"]["coordinates"],
            dtype=float
        )
        coords = coords.reshape(-1, coords.shape[-1])[:, :2]

        moved = np.min(
            np.hypot(
                (coords[:, 0] - lon) * np.cos(np.radians(lat)),
                coords[:, 1] - lat
            )
        ) * 111.0

        print(
            "      seed (%.4f, %.4f) -> comid %s, snapped %.2f km"
            % (lat, lon, comid, moved)
        )

    except Exception:
        print(
            "      seed (%.4f, %.4f) -> comid %s"
            % (lat, lon, comid)
        )

    return comid


def seed_path(seed):
    if seed.get("site"):
        return (
            "nwissite/USGS-%s" % seed["site"],
            "site %s" % seed["site"]
        )

    if seed.get("comid"):
        return (
            "comid/%s" % seed["comid"],
            "comid %s" % seed["comid"]
        )

    if seed.get("point"):
        comid = comid_at_point(*seed["point"])

        if comid is None:
            return None, None

        return (
            "comid/%s" % comid,
            "comid %s (from point)" % comid
        )

    return None, None


def fetch_flowlines_seeded(name, prefix, mode, distance):
    cached = read_cached_geojson(name)

    if cached is not None:
        print("   %-22s cached (%d features)" % (name, len(cached)))
        return cached

    if prefix is None:
        return None

    data = nldi_get(
        "%s/navigation/%s/flowlines" % (prefix, mode),
        {"f": "json", "distance": distance}
    )

    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None

    frame = gpd.GeoDataFrame.from_features(data["features"], crs=WGS84)

    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")

    print("   %-22s fetched (%d features)" % (name, len(frame)))

    return frame


def fetch_local_tributary(trib, mainstem_comids):
    label = trib["label"]
    stem = label.split()[0].lower()

    up_name = "local_%s_um.geojson" % stem
    down_name = "local_%s_dm.geojson" % stem

    up = read_cached_geojson(up_name)
    down = read_cached_geojson(down_name)

    note = "cached"

    if up is None or down is None:
        prefix, note = seed_path(trib["seed"])

        if prefix is None:
            print("   %s NOT DRAWN (seed did not resolve)" % label)
            return None

        up = fetch_flowlines_seeded(
            up_name,
            prefix,
            "UM",
            trib.get("up_km", 40)
        )

        down = fetch_flowlines_seeded(
            down_name,
            prefix,
            "DM",
            trib.get("down_km", 15)
        )

    down = drop_comids(down, mainstem_comids)

    parts = [
        f for f in (up, down)
        if f is not None and len(f)
    ]

    if not parts:
        print("   %s NOT DRAWN (no flowlines)" % label)
        return None

    frame = gpd.GeoDataFrame(
        pd.concat(parts, ignore_index=True),
        crs=WGS84
    )

    print(
        "   %-22s %d features (seed %s)"
        % (label, len(frame), note)
    )

    return frame


def fetch_basin_seeded(name, seed, fallback_site):
    cached = read_cached_geojson(name)

    if cached is None:
        prefix, note = seed_path(seed)

        data = (
            nldi_get("%s/basin" % prefix, {"f": "json"})
            if prefix else None
        )

        if not data:
            print(
                "   basin from seed UNAVAILABLE -- falling back to the "
                "basin above gage %s" % fallback_site
            )
            return fetch_basin(name, fallback_site)

        cached = gpd.GeoDataFrame.from_features(
            data["features"],
            crs=WGS84
        )

        ensure_cache_dir()
        cached.to_file(cache_path(name), driver="GeoJSON")

        print(
            "   %-22s fetched (seed %s)"
            % (name, note)
        )

    else:
        print("   %-22s cached" % name)

    try:
        area = (
            cached.to_crs("EPSG:5070").area.sum()
            / 2.589988e6
        )

        flag = (
            ""
            if abs(area - BASIN_EXPECT_SQ_MI) / BASIN_EXPECT_SQ_MI < 0.15
            else "   <-- CHECK: expected about %.0f" % BASIN_EXPECT_SQ_MI
        )

        print(
            "   %-22s %.0f sq mi%s"
            % ("basin area", area, flag)
        )

    except Exception:
        pass

    return cached


def fetch_basin(name, site):
    cached = read_cached_geojson(name)

    if cached is not None:
        print("   %-22s cached" % name)
        return cached

    data = nldi_get(
        "nwissite/USGS-%s/basin" % site,
        {"f": "json"}
    )

    if not data:
        print("   %-22s UNAVAILABLE" % name)
        return None

    frame = gpd.GeoDataFrame.from_features(
        data["features"],
        crs=WGS84
    )

    ensure_cache_dir()
    frame.to_file(cache_path(name), driver="GeoJSON")

    print("   %-22s fetched" % name)

    return frame


def fetch_sites(numbers):
    path = cache_path("sites.csv")

    if os.path.isfile(path):
        table = pd.read_csv(path, dtype={"site_no": str})

        has_fallback = bool(
            table.get(
                "is_fallback",
                pd.Series([False])
            ).any()
        )

        if set(numbers) <= set(table["site_no"]) and not has_fallback:
            print("   %-22s cached" % "gage coordinates")
            return table, False

        if has_fallback:
            print(
                "   %-22s cached copy holds fallbacks -- refetching"
                % "gage coordinates"
            )

    print("   NWIS %s" % NWIS_SITE_URL)

    text = http_get(
        NWIS_SITE_URL,
        {
            "format": "rdb",
            "sites": ",".join(numbers),
            "siteOutput": "expanded"
        }
    )

    rows = []

    if text:
        lines = [
            ln for ln in text.splitlines()
            if ln and not ln.startswith("#")
        ]

        if len(lines) >= 3:
            header = lines[0].split("\t")

            for line in lines[2:]:
                parts = line.split("\t")

                if len(parts) != len(header):
                    continue

                rec = dict(zip(header, parts))

                try:
                    rows.append(
                        {
                            "site_no": rec["site_no"].strip(),
                            "station_nm": rec.get(
                                "station_nm", ""
                            ).strip(),
                            "lat": float(rec["dec_lat_va"]),
                            "lon": float(rec["dec_long_va"]),
                            "is_fallback": False
                        }
                    )

                except (KeyError, ValueError):
                    continue

    got = {r["site_no"] for r in rows}
    missing = [n for n in numbers if n not in got]

    for number in missing:
        if number in FALLBACK_SITES:
            lat, lon = FALLBACK_SITES[number]

            rows.append(
                {
                    "site_no": number,
                    "station_nm": "(fallback)",
                    "lat": lat,
                    "lon": lon,
                    "is_fallback": True
                }
            )

    if missing:
        print(
            "   *** NWIS did not return %s -- using FALLBACK_SITES, which are"
            % ", ".join(missing)
        )
        print(
            "       approximate. The figure is stamped accordingly."
        )

    table = pd.DataFrame(rows)

    if len(table):
        ensure_cache_dir()
        table.to_csv(path, index=False)

    print(
        "   %-22s %d site(s)"
        % ("gage coordinates", len(table))
    )

    return (
        table,
        bool(len(table))
        and bool(table["is_fallback"].any())
    )


def to_web(frame):
    return (
        None
        if frame is None or not len(frame)
        else frame.to_crs(WEBM)
    )


def clip_to(frame, basin):
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
        geometry=[
            Point(x, y)
            for x, y in zip(table["lon"], table["lat"])
        ],
        crs=WGS84
    ).to_crs(WEBM)


def longest_strand(frame):
    from shapely.ops import linemerge

    if frame is None or not len(frame):
        return None

    merged = linemerge(list(frame.geometry))

    line = (
        max(merged.geoms, key=lambda g: g.length)
        if merged.geom_type == "MultiLineString"
        else merged
    )

    return line if line.length > 0 else None


def label_along_line(
    ax,
    frame,
    text,
    frac,
    color,
    size=10,
    offset=(0, 0),
    rotation_override=None
):
    line = longest_strand(frame)

    if line is None:
        return

    point = line.interpolate(frac, normalized=True)

    rotation = 0.0

    if ROTATE_RIVER_LABELS:
        step = 0.04

        before = line.interpolate(
            max(frac - step, 0.0),
            normalized=True
        )

        after = line.interpolate(
            min(frac + step, 1.0),
            normalized=True
        )

        rotation = np.degrees(
            np.arctan2(
                after.y - before.y,
                after.x - before.x
            )
            
        )

        if rotation > 90:
            rotation -= 180
        elif rotation < -90:
            rotation += 180
        
    if rotation_override is not None:
        rotation = rotation_override

    ax.annotate(
        text,
        xy=(point.x, point.y),
        xytext=offset,
        textcoords="offset points",
        rotation=rotation,
        rotation_mode="anchor",
        ha="center",
        va="bottom",
        zorder=11,
        style="italic",
        **halo(size, "bold", color)
    )


def annotate_leader(ax, xy, text, offset, color, size=9):
    far = (
        offset[0] ** 2
        + offset[1] ** 2
    ) ** 0.5 >= LEADER_MIN_POINTS

    ax.annotate(
        text,
        xy=xy,
        xytext=offset,
        textcoords="offset points",
        ha="left" if offset[0] >= 0 else "right",
        va="bottom" if offset[1] >= 0 else "top",
        zorder=11,
        arrowprops=(
            dict(
                arrowstyle="-",
                color=color,
                lw=0.9,
                shrinkA=0,
                shrinkB=6,
                alpha=0.85
            )
            if far else None
        ),
        **halo(size, "bold", color)
    )


def halo(size, weight="normal", color="black"):
    import matplotlib.patheffects as pe

    return dict(
        fontsize=size,
        fontweight=weight,
        color=color,
        path_effects=[
            pe.withStroke(
                linewidth=LABEL_HALO,
                foreground="white"
            )
        ]
    )


def add_scalebar(ax, length_mi=None):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    length_mi = (
        SCALEBAR_MI
        if length_mi is None
        else length_mi
    )

    lat = np.degrees(
        2 * np.arctan(
            np.exp(
                (y0 + y1)
                / 2
                / 6378137.0
            )
        )
        - np.pi / 2
    )

    metres = (
        length_mi
        * 1609.34
        / np.cos(np.radians(lat))
    )

    x_start = (
        x0
        + SCALEBAR_ANCHOR[0]
        * (x1 - x0)
    )

    y_pos = (
        y0
        + SCALEBAR_ANCHOR[1]
        * (y1 - y0)
    )

    ax.plot(
        [x_start, x_start + metres],
        [y_pos, y_pos],
        color="black",
        lw=3,
        solid_capstyle="butt",
        zorder=12
    )

    for x in (x_start, x_start + metres):
        ax.plot(
            [x, x],
            [
                y_pos,
                y_pos + 0.012 * (y1 - y0)
            ],
            color="black",
            lw=3,
            zorder=12
        )

    ax.text(
        x_start + metres / 2,
        y_pos + 0.018 * (y1 - y0),
        "%d mi" % length_mi,
        ha="center",
        va="bottom",
        **halo(9, "bold")
    )


def add_north_arrow(ax):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    x = x0 + 0.045 * (x1 - x0)
    y = y1 - 0.135 * (y1 - y0)

    ax.annotate(
        "N",
        xy=(x, y),
        xytext=(
            x,
            y - 0.055 * (y1 - y0)
        ),
        ha="center",
        va="center",
        zorder=12,
        arrowprops=dict(
            arrowstyle="-|>",
            color="black",
            lw=2
        ),
        **halo(12, "bold")
    )


def main():
    ensure_cache_dir()

    print("=" * 78)
    print(
        "Basemap   : Esri.%s (zoom %s)"
        % (BASEMAP, ZOOM)
    )
    print(
        "Cache     : %s"
        % os.path.abspath(CACHE_DIR)
    )
    print("=" * 78)

    sites, used_fallback = fetch_sites(
        [n for n, _, _, _ in GAGES]
    )

    basin = fetch_basin_seeded(
        "basin.geojson",
        BASIN_SEED,
        OUTLET_SITE
    )

    network = fetch_flowlines(
        "network_ut.geojson",
        OUTLET_SITE,
        "UT"
    )

    mainstem = fetch_flowlines(
        "cowlitz_um.geojson",
        OUTLET_SITE,
        "UM"
    )

    mainstem_comids = comid_set(mainstem)

    if mainstem is not None and not mainstem_comids:
        print(
            "   NOTE: no comid column on the mainstem, so the tributaries "
            "cannot be"
        )
        print(
            "         trimmed at their confluences and will run on down the "
            "Cowlitz."
        )

    tribs = [
        (
            t,
            fetch_tributary(
                t,
                mainstem_comids
            )
        )
        for t in TRIBUTARIES
    ]

    downstream = fetch_flowlines_seeded(
        "cowlitz_dm.geojson",
        "nwissite/USGS-%s" % OUTLET_SITE,
        "DM",
        DOWNSTREAM_KM
    )

    local_trim = (
        mainstem_comids
        | comid_set(downstream)
    )

    locals_ = [
        (
            t,
            fetch_local_tributary(
                t,
                local_trim
            )
        )
        for t in LOCAL_TRIBUTARIES
    ]

    if all(
        x is None
        for x in (
            [basin, network, mainstem, downstream]
            + [f for _, f in tribs]
            + [f for _, f in locals_]
        )
    ):
        raise SystemExit(
            "No geometry was fetched and the cache is empty, so there is "
            "nothing to draw.\nCheck the network, then re-run. The gage "
            "coordinates alone are not a map."
        )

    basin_w = to_web(basin)
    network_w = clip_to(to_web(network), basin_w)
    mainstem_w = clip_to(to_web(mainstem), basin_w)

    tribs_w = [
        (
            t,
            clip_to(to_web(f), basin_w)
        )
        for t, f in tribs
    ]

    downstream_w = to_web(downstream)

    locals_w = [
        (
            t,
            to_web(f)
        )
        for t, f in locals_
    ]

    fig, ax = plt.subplots(
        figsize=FIGSIZE
    )

    extent_parts = [
        f
        for f in (
            [basin_w, network_w, downstream_w]
            + [f for _, f in locals_w]
        )
        if f is not None and len(f)
    ]

    if not extent_parts:
        raise SystemExit(
            "nothing to set an extent from"
        )

    bounds = np.array(
        [f.total_bounds for f in extent_parts]
    )

    minx = bounds[:, 0].min()
    miny = bounds[:, 1].min()
    maxx = bounds[:, 2].max()
    maxy = bounds[:, 3].max()

    dx = (
        (maxx - minx)
        * EXTENT_PAD
    )

    dy = (
        (maxy - miny)
        * EXTENT_PAD
    )

    ax.set_xlim(
        minx - dx,
        maxx + dx
    )

    ax.set_ylim(
        miny - dy,
        maxy + dy
    )

    if basin_w is not None:
        basin_w.plot(
            ax=ax,
            facecolor=C_BASIN,
            alpha=BASIN_FILL_ALPHA,
            edgecolor="none",
            zorder=2
        )

        basin_w.boundary.plot(
            ax=ax,
            color=C_BASIN,
            lw=1.6,
            ls="--",
            zorder=3
        )

    if network_w is not None:
        network_w.plot(
            ax=ax,
            color=C_NET,
            lw=LW_NET,
            zorder=4
        )

    if mainstem_w is not None:
        mainstem_w.plot(
            ax=ax,
            color=C_MAIN,
            lw=LW_MAIN,
            zorder=6
        )

    if downstream_w is not None and len(downstream_w):
        downstream_w.plot(
            ax=ax,
            color=C_MAIN,
            lw=LW_MAIN,
            zorder=6
        )

    for _, frame in locals_w:
        if frame is not None and len(frame):
            frame.plot(
                ax=ax,
                color=C_TRIB,
                lw=LW_TRIB,
                zorder=5
            )

    for _, frame in tribs_w:
        if frame is not None and len(frame):
            frame.plot(
                ax=ax,
                color=C_TRIB,
                lw=LW_TRIB,
                zorder=5
            )

    # --- river names ----------------------------------------------------------
    for trib, frame in tribs_w:
        label_along_line(
            ax,
            frame,
            trib["label"],
            trib.get(
                "label_frac",
                0.45
            ),
            C_TRIB
        )

    for trib, frame in locals_w:
        label_along_line(
            ax,
            frame,
            trib["label"],
            trib.get(
                "label_frac",
                0.45
            ),
            C_TRIB,
            offset=trib.get(
                "label_offset",
                (0, 0)
            ),
            rotation_override=trib.get(
                "label_rotation"
            )
        )

    if (
        MAINSTEM_LABEL
        and mainstem_w is not None
    ):
        label_along_line(
            ax,
            mainstem_w,
            MAINSTEM_LABEL,
            MAINSTEM_LABEL_FRAC,
            C_MAIN,
            size=11
        )

    # --- gages ---------------------------------------------------------------
    gage_pts = (
        points_frame(
            sites[
                sites["site_no"].isin(
                    [
                        n
                        for n, _, _, _
                        in GAGES
                    ]
                )
            ]
        )
        if len(sites)
        else None
    )

    if gage_pts is not None:
        by_site = {
            r["site_no"]: r
            for _, r
            in gage_pts.iterrows()
        }

        for (
            number,
            label,
            sub,
            offset
        ) in GAGES:

            if number not in by_site:
                print(
                    "   gage %s has no coordinate -- not drawn"
                    % number
                )
                continue

            row = by_site[number]

            ax.plot(
                [row.geometry.x],
                [row.geometry.y],
                marker="o",
                ms=11,
                mfc=C_GAGE,
                mec="white",
                mew=1.8,
                zorder=10
            )

            annotate_leader(
                ax,
                (
                    row.geometry.x,
                    row.geometry.y
                ),
                "%s\n%s"
                % (label, sub),
                offset,
                C_GAGE
            )

    # --- dams ----------------------------------------------------------------
    if SHOW_DAMS:
        dam_pts = points_frame(
            pd.DataFrame(
                [
                    {
                        "lat": la,
                        "lon": lo
                    }
                    for _, la, lo, _, _
                    in DAMS
                ]
            )
        )

        for (
            name,
            _,
            _,
            note,
            offset
        ), (_, row) in zip(
            DAMS,
            dam_pts.iterrows()
        ):

            ax.plot(
                [row.geometry.x],
                [row.geometry.y],
                marker="s",
                ms=10,
                mfc=C_DAM,
                mec="white",
                mew=1.6,
                zorder=10
            )

            annotate_leader(
                ax,
                (
                    row.geometry.x,
                    row.geometry.y
                ),
                "%s\n%s"
                % (name, note),
                offset,
                C_DAM
            )

    # --- lake labels ---------------------------------------------------------
    if SHOW_LAKE_LABELS:
        lake_pts = points_frame(
            pd.DataFrame(
                [
                    {
                        "lat": la,
                        "lon": lo
                    }
                    for la, lo, _, _
                    in LAKES
                ]
            )
        )

        for (
            _,
            _,
            label,
            rot
        ), (_, row) in zip(
            LAKES,
            lake_pts.iterrows()
        ):

            ax.text(
                row.geometry.x,
                row.geometry.y,
                label,
                rotation=rot,
                ha="center",
                va="center",
                style="italic",
                zorder=11,
                **halo(
                    9.5,
                    "normal",
                    "#12507a"
                )
            )

    # --- basemap -------------------------------------------------------------
    try:
        if not hasattr(
            cx.providers.Esri,
            BASEMAP
        ):
            raise SystemExit(
                "BASEMAP is '%s', which this xyzservices does not have.\n"
                "Available Esri basemaps: %s"
                % (
                    BASEMAP,
                    ", ".join(
                        sorted(
                            cx.providers.Esri.keys()
                        )
                    )
                )
            )

        source = getattr(
            cx.providers.Esri,
            BASEMAP
        )

        cx.add_basemap(
            ax,
            source=source,
            zoom=ZOOM,
            alpha=BASEMAP_ALPHA,
            attribution_size=6
        )

        if (
            BASEMAP == "WorldImagery"
            and IMAGERY_LABEL_OVERLAY
        ):
            cx.add_basemap(
                ax,
                source=ESRI_REFERENCE_URL,
                zoom=ZOOM,
                attribution=False
            )

    except SystemExit:
        raise

    except Exception as exc:
        print(
            "\n   *** BASEMAP NOT DRAWN: %s"
            % exc
        )

        print(
            "   The vectors are still correct; only the backdrop is missing."
        )

        print(
            "   Usually the tile host is unreachable. contextily caches "
            "tiles, so a"
        )

        print(
            "   later run on a working connection will fill it in."
        )

        ax.set_facecolor(
            "#eef2f5"
        )

    ax.set_xlim(
        minx - dx,
        maxx + dx
    )

    ax.set_ylim(
        miny - dy,
        maxy + dy
    )

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_edgecolor(
            "#33414d"
        )
        spine.set_linewidth(
            1.0
        )

    if SHOW_SCALEBAR:
        add_scalebar(ax)

    if SHOW_NORTH_ARROW:
        add_north_arrow(ax)

    ax.set_title(
        TITLE,
        fontsize=13.5,
        fontweight="bold",
        pad=10
    )

    note = SUBTITLE

    if used_fallback:
        note += (
            "\nWARNING: one or more gage positions came from "
            "FALLBACK_SITES (approximate), not from NWIS."
        )

    ax.text(
        0.5,
        -0.035,
        note,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.5,
        color=(
            "#c0392b"
            if used_fallback
            else "#4a5568"
        )
    )

    fig.savefig(
        OUT_PNG,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(fig)

    print(
        "\nWrote %s  (%.1f x %.1f in at %d dpi)"
        % (
            os.path.abspath(OUT_PNG),
            FIGSIZE[0],
            FIGSIZE[1],
            DPI
        )
    )

    if used_fallback:
        print(
            "*** At least one gage used FALLBACK_SITES. Delete %s and re-run"
            % cache_path("sites.csv")
        )

        print(
            "    with a working connection to get the published coordinates."
        )


main()