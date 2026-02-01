#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
obfuscate_bathy_tiles.py

Goal
----
Generate (optionally obfuscated) bathymetric contour vector tiles (MVT) packed into MBTiles,
starting from EMODnet Bathymetry DTM tiles.

What this script does
---------------------
1) Given a geographic bounding box, it discovers which EMODnet DTM tiles intersect it (via CSW metadata).
2) Downloads only those tiles (ZIP), extracts GeoTIFF(s), and caches them locally.
3) For each zoom:
   - Builds a mosaic of the input tiles (in the source CRS),
   - Warps + clips it into WebMercator (EPSG:3857) at a zoom-appropriate resolution,
   - Writes that as a cached intermediate raster GeoTIFF (fast for repeated tile reads).
4) For each WebMercator tile (z/x/y) inside the bbox:
   - Reads a tile-sized raster window from the intermediate raster,
   - Builds a land mask (OSM land polygons) and drops contours on land,
   - Generates bathymetric contours at the configured interval for that zoom,
   - Applies optional obfuscation (jitter/simplify/drop/level-noise),
   - Encodes vector features into an MVT tile and writes to MBTiles.
5) Optionally stores per-tile bathymetry “payload” for client-side processing using one of several modes:
   - none
   - mvt_grid_i16  (quantized int16 grid, zlib+base64 in MVT properties)
   - mvt_grid_f32  (float32 grid, zlib+base64 in MVT properties)
   - mvt_png       (small grayscale PNG, zlib+base64 in MVT properties)

Install
-------
pip install -r requirements.txt

Run
---
python obfuscate_bathy_tiles.py config.json
"""

# ----------------------------
# Standard library imports
# ----------------------------
from __future__ import annotations  # Python 3.7+ typing improvements

import base64  # For embedding binary payloads in MVT properties as ASCII
import json  # Configuration file parsing
import math  # Numeric utilities
import os  # Path manipulation
import random  # Deterministic obfuscation with seed
import re  # Lightweight parsing of CSW XML for download URLs
import shutil  # Move extracted files, cleanup
import sqlite3  # MBTiles container is an SQLite DB
import sys  # CLI args, exit codes
import time  # Simple timing/logging
import zipfile  # EMODnet DTM tiles typically distributed as .zip
import zlib  # Compress grid payloads (fast and widely supported)

from dataclasses import dataclass  # Clean config model
from typing import Any, Dict, Iterable, List, Optional, Tuple  # Type hints

# ----------------------------
# Third-party imports
# ----------------------------
import mercantile  # Slippy map tile math (z/x/y) + bounds in lon/lat
import numpy as np  # Fast array operations
import requests  # HTTP requests for CSW and tile downloads

import rasterio  # Raster I/O
import rasterio.warp  # Reprojection utilities
from rasterio.enums import Resampling  # Resampling methods
from rasterio.features import rasterize  # Rasterize polygons into masks
from rasterio.merge import merge as rio_merge  # Mosaic multiple rasters
from rasterio.transform import from_bounds  # Build transforms from bbox+size

import fiona  # Read shapefiles (OSM land polygons)
from shapely.geometry import LineString, Polygon, MultiPolygon, box, shape  # Geometry objects
from shapely.ops import clip_by_rect  # Fast rectangular clipping
from shapely.strtree import STRtree  # Spatial index for land polygons

# Matplotlib is used only for contour extraction (robust and widely available)
import matplotlib
matplotlib.use("Agg")  # Headless backend (no GUI needed)
import matplotlib.pyplot as plt

import mapbox_vector_tile  # MVT encoding


# ============================================================
# 1) Small helpers: zoom-range configuration
# ============================================================

def _parse_zoom_key(key: str) -> Tuple[int, int]:
    """
    Parse zoom key strings like:
      "12"    -> (12, 12)
      "10-14" -> (10, 14)
    """
    key = key.strip()  # Defensive whitespace removal
    if "-" in key:  # Range case
        a, b = key.split("-", 1)
        return int(a), int(b)
    z = int(key)  # Single zoom
    return z, z


def _value_for_zoom(ranged: Dict[str, Any], z: int, default: Any = None) -> Any:
    """
    Get a value for the specific zoom `z` from a dict whose keys can be:
      - exact zoom ("12")
      - zoom range ("10-14")
    Exact zoom entries take priority over ranges.
    """
    exact = str(z)  # Exact key string
    if exact in ranged:
        return ranged[exact]  # Exact match wins

    # Otherwise, scan ranges and pick first that matches
    for k, v in ranged.items():
        z0, z1 = _parse_zoom_key(k)
        if z0 <= z <= z1:
            return v
    return default  # Nothing matched


# ============================================================
# 2) MBTiles writer
# ============================================================

# Minimal MBTiles schema for vector tiles:
MBTILES_SCHEMA = """
CREATE TABLE IF NOT EXISTS metadata (name TEXT, value TEXT);
CREATE TABLE IF NOT EXISTS tiles (zoom_level INTEGER, tile_column INTEGER, tile_row INTEGER, tile_data BLOB);
CREATE UNIQUE INDEX IF NOT EXISTS tile_index ON tiles (zoom_level, tile_column, tile_row);
"""

def mbtiles_open(path: str) -> sqlite3.Connection:
    """
    Create/open an MBTiles file (SQLite) and ensure required schema exists.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)  # Ensure output folder exists
    conn = sqlite3.connect(path)  # Create/open DB
    conn.execute("PRAGMA journal_mode=WAL;")  # WAL improves write performance/safety
    conn.executescript(MBTILES_SCHEMA)  # Create tables if missing
    return conn


def mbtiles_put_metadata(conn: sqlite3.Connection, md: Dict[str, str]) -> None:
    """
    Overwrite MBTiles metadata table.
    """
    conn.execute("DELETE FROM metadata;")  # Reset metadata
    for k, v in md.items():
        conn.execute("INSERT INTO metadata(name, value) VALUES(?, ?);", (k, v))
    conn.commit()  # Persist


def mbtiles_put_tile(conn: sqlite3.Connection, z: int, x: int, y_tms: int, data: bytes) -> None:
    """
    Write one tile blob into MBTiles.
    Note: MBTiles tile_row uses TMS (flipped Y).
    """
    conn.execute(
        "INSERT OR REPLACE INTO tiles(zoom_level, tile_column, tile_row, tile_data) VALUES(?,?,?,?);",
        (z, x, y_tms, data),
    )


def mbtiles_commit(conn: sqlite3.Connection) -> None:
    """
    Commit pending inserts.
    """
    conn.commit()


def tms_y(z: int, y_xyz: int) -> int:
    """
    Convert XYZ y to TMS y (MBTiles uses TMS row indexing).
    """
    return (2**z - 1) - y_xyz


# ============================================================
# 3) WebMercator math
# ============================================================

WEBMERC_MAX = 20037508.342789244  # Half-world extent in EPSG:3857 meters

def lonlat_bbox_to_merc(b: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    """
    Transform lon/lat bbox (EPSG:4326) to WebMercator bbox (EPSG:3857).
    """
    minlon, minlat, maxlon, maxlat = b
    tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)  # Always treat inputs as (lon,lat)
    minx, miny = tf.transform(minlon, minlat)
    maxx, maxy = tf.transform(maxlon, maxlat)
    return (minx, miny, maxx, maxy)


def tile_bounds_merc(z: int, x: int, y: int) -> Tuple[float, float, float, float]:
    """
    Get WebMercator bounds of a slippy tile z/x/y.
    """
    b = mercantile.bounds(x, y, z)  # Returns lon/lat bounds
    return lonlat_bbox_to_merc((b.west, b.south, b.east, b.north))


def meters_per_pixel_equator(z: int, tile_size: int = 256) -> float:
    """
    Approximate resolution (meters/pixel) at the equator for WebMercator tiles.
    """
    world = 2 * WEBMERC_MAX  # Total world width in meters
    return world / (tile_size * 2**z)  # Meters per pixel at equator


def zoom_target_resolution_m(z: int, tile_pixels: int, tile_size: int = 256) -> float:
    """
    If we later read tiles into tile_pixels x tile_pixels, pick a source resolution
    that roughly matches that target.
    """
    # If the canonical tile is 256px but we read 512px, we want half the m/px.
    return meters_per_pixel_equator(z, tile_size=tile_size) * (tile_size / tile_pixels)


# ============================================================
# 4) EMODnet tile discovery via CSW (GeoNetwork)
# ============================================================

# CSW query template for bounding-box search of EMODnet DTM tile records.
CSW_GETRECORDS_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<csw:GetRecords
  xmlns:csw="http://www.opengis.net/cat/csw/2.0.2"
  xmlns:ogc="http://www.opengis.net/ogc"
  service="CSW"
  version="2.0.2"
  resultType="results"
  startPosition="{start}"
  maxRecords="{max_records}"
  outputSchema="http://www.isotc211.org/2005/gmd"
  outputFormat="application/xml">
  <csw:Query typeNames="csw:Record">
    <csw:ElementSetName>full</csw:ElementSetName>
    <csw:Constraint version="1.1.0">
      <ogc:Filter>
        <ogc:And>
          <ogc:PropertyIsLike wildCard="*" singleChar="?" escapeChar="\\">
            <ogc:PropertyName>AnyText</ogc:PropertyName>
            <ogc:Literal>*EMODnet Digital Bathymetry (DTM {year}) - Tile*</ogc:Literal>
          </ogc:PropertyIsLike>
          <ogc:BBOX>
            <ogc:PropertyName>ows:BoundingBox</ogc:PropertyName>
            <gml:Envelope xmlns:gml="http://www.opengis.net/gml">
              <gml:lowerCorner>{minx} {miny}</gml:lowerCorner>
              <gml:upperCorner>{maxx} {maxy}</gml:upperCorner>
            </gml:Envelope>
          </ogc:BBOX>
        </ogc:And>
      </ogc:Filter>
    </csw:Constraint>
  </csw:Query>
</csw:GetRecords>
"""

def csw_search_tiles(
    csw_url: str,
    bbox_lonlat: Tuple[float, float, float, float],
    year: int,
    timeout_s: int = 60,
    max_pages: int = 50,
) -> List[str]:
    """
    Query EMODnet GeoNetwork CSW for DTM tile records intersecting bbox.
    Then extract download URLs from returned XML with a regex.

    Returns a list of URLs like:
      https://downloads.emodnet-bathymetry.eu/.../TILE_YYYY.tif.zip

    Notes:
    - CSW parsing is intentionally lightweight (regex) to minimize dependencies.
    - If EMODnet metadata structure changes, you may need to adjust the regex/query.
    """
    minx, miny, maxx, maxy = bbox_lonlat  # CSW bbox is lon/lat

    tile_urls: List[str] = []  # Output list
    seen = set()  # De-duplication

    # Regex to find URLs ending with _{year}.tif.zip under downloads.emodnet-bathymetry.eu
    url_re = re.compile(r"https://downloads\.emodnet-bathymetry\.eu/[^\s\"<]+?_{}\.tif\.zip".format(year))

    start = 1  # CSW paging: startPosition is 1-based
    max_records = 50  # Page size
    headers = {"Content-Type": "application/xml"}  # CSW POST body is XML

    for _ in range(max_pages):  # Safety limit on paging
        body = CSW_GETRECORDS_TEMPLATE.format(
            start=start,
            max_records=max_records,
            year=year,
            minx=minx,
            miny=miny,
            maxx=maxx,
            maxy=maxy,
        )

        # POST the CSW request
        r = requests.post(csw_url, data=body.encode("utf-8"), headers=headers, timeout=timeout_s)
        r.raise_for_status()  # Throw if HTTP error
        xml = r.text  # Raw XML as text

        # Extract URLs matching our regex
        found = url_re.findall(xml)
        for u in found:
            if u not in seen:
                seen.add(u)
                tile_urls.append(u)

        # Attempt to stop when no more records are available
        if "nextRecord=\"0\"" in xml or "numberOfRecordsReturned=\"0\"" in xml:
            break
        if len(found) == 0 and start > 1:
            break

        start += max_records  # Next page

    return tile_urls


# ============================================================
# 5) Download + unzip cache
# ============================================================

def download_file(url: str, out_path: str, timeout_s: int = 180) -> None:
    """
    Download a URL to a local file path (streamed).
    If file exists, it is not downloaded again.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)  # Ensure directory exists
    if os.path.exists(out_path):
        return  # Cached already

    # Stream download to avoid huge memory usage
    with requests.get(url, stream=True, timeout=timeout_s) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):  # 1MB chunks
                if chunk:
                    f.write(chunk)


def download_and_extract_zip(url: str, out_dir: str, timeout_s: int = 180) -> str:
    """
    Download a .zip and extract the first GeoTIFF inside.

    Returns: path to extracted GeoTIFF.
    """
    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.basename(url)  # Local zip file name
    local_zip = os.path.join(out_dir, fname)

    download_file(url, local_zip, timeout_s=timeout_s)  # Download if needed

    # Open the zip and find a .tif/.tiff member
    with zipfile.ZipFile(local_zip, "r") as z:
        tif_members = [m for m in z.namelist() if m.lower().endswith((".tif", ".tiff"))]
        if not tif_members:
            raise RuntimeError(f"No GeoTIFF found inside {local_zip}")

        member = tif_members[0]  # Use first GeoTIFF in the archive
        out_tif = os.path.join(out_dir, os.path.basename(member))  # Flatten into out_dir

        if not os.path.exists(out_tif):
            z.extract(member, out_dir)  # Extract possibly into subfolders
            extracted = os.path.join(out_dir, member)

            # If it extracted under a subfolder, move it to out_tif
            if extracted != out_tif:
                os.makedirs(os.path.dirname(out_tif), exist_ok=True)
                shutil.move(extracted, out_tif)

        return out_tif


# ============================================================
# 6) Coastline/land mask via OSM land polygons
# ============================================================

def ensure_osm_land_polygons(cache_dir: str, url: str) -> str:
    """
    Download and extract OSM land-polygons zip into cache, if needed.

    Returns: root directory where shapefiles live.
    """
    coast_dir = os.path.join(cache_dir, "coastline_mask", "osm_land_polygons")  # Cache location
    os.makedirs(coast_dir, exist_ok=True)

    zip_path = os.path.join(coast_dir, os.path.basename(url))  # Local zip
    download_file(url, zip_path, timeout_s=600)  # Land polygons can be large; allow more time

    marker = os.path.join(coast_dir, ".extracted")  # Marker file to avoid re-extract
    if not os.path.exists(marker):
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(coast_dir)  # Extract entire dataset
        with open(marker, "w", encoding="utf-8") as f:
            f.write("ok\n")

    return coast_dir


def find_shapefiles(root_dir: str) -> List[str]:
    """
    Recursively find all .shp files under root_dir.
    """
    shps: List[str] = []
    for base, _, files in os.walk(root_dir):
        for fn in files:
            if fn.lower().endswith(".shp"):
                shps.append(os.path.join(base, fn))
    return shps


def bbox_intersects(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> bool:
    """
    Axis-aligned bbox intersection test.
    """
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])


def load_land_polygons_strtree(
    land_polygons_root: str,
    bbox_lonlat: Tuple[float, float, float, float],
    max_features: Optional[int] = None,
) -> Tuple[STRtree, List[Polygon | MultiPolygon]]:
    """
    Load land polygons overlapping bbox_lonlat (EPSG:4326) and build an STRtree index.

    Returns:
      (tree, geometries)

    Notes:
    - The "land-polygons-split-4326" dataset is in EPSG:4326, which matches bbox_lonlat.
    - For huge bboxes this can be heavy; you can cap `max_features` for testing.
    """
    shps = find_shapefiles(land_polygons_root)
    if not shps:
        raise RuntimeError(f"No shapefiles found under {land_polygons_root}")

    geoms: List[Polygon | MultiPolygon] = []
    count = 0

    for shp in shps:
        with fiona.open(shp, "r") as src:
            # shapefile bounds are in its CRS (here: EPSG:4326)
            b = src.bounds  # (minx, miny, maxx, maxy)
            if not bbox_intersects(b, bbox_lonlat):
                continue  # Skip shapefiles fully outside bbox

            # Read features and keep polygonal geometries
            for feat in src:
                g = feat.get("geometry")
                if not g:
                    continue
                s = shape(g)  # Convert GeoJSON-like mapping to shapely geometry
                if s.is_empty:
                    continue
                if isinstance(s, (Polygon, MultiPolygon)):
                    geoms.append(s)
                    count += 1
                    if max_features and count >= max_features:
                        break

        if max_features and count >= max_features:
            break

    # If bbox is offshore you may have zero land features; that's OK.
    tree = STRtree(geoms) if geoms else STRtree([])
    return tree, geoms


def rasterize_land_mask(
    land_tree: STRtree,
    land_geoms: List[Polygon | MultiPolygon],
    tile_bounds_merc: Tuple[float, float, float, float],
    out_shape: Tuple[int, int],
) -> np.ndarray:
    """
    Rasterize land polygons into a boolean mask in tile pixel grid.

    Input:
      - land_geoms in EPSG:4326
      - tile_bounds_merc in EPSG:3857
      - out_shape = (H, W) of the raster we will contour

    Output:
      - mask[H,W] True where land, False elsewhere.

    Method:
      1) Convert tile bounds to lon/lat bbox
      2) Query STRtree in lon/lat space
      3) Reproject candidates to EPSG:3857
      4) Rasterize onto the tile pixel grid
    """
    if not land_geoms:
        return np.zeros(out_shape, dtype=bool)  # No land => all water

    # Transform tile bounds from EPSG:3857 to EPSG:4326 to query the STRtree
    tf_back = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    minx, miny, maxx, maxy = tile_bounds_merc
    w, s = tf_back.transform(minx, miny)
    e, n = tf_back.transform(maxx, maxy)

    tile_bbox_lonlat = (min(w, e), min(s, n), max(w, e), max(s, n))  # Ensure proper min/max
    tile_box_lonlat = box(*tile_bbox_lonlat)  # Shapely bbox polygon

    # STRtree query: returns geometries whose bbox intersects tile_box_lonlat
    candidates = land_tree.query(tile_box_lonlat)
    if not candidates:
        return np.zeros(out_shape, dtype=bool)  # No land nearby

    # Transformer forward: EPSG:4326 -> EPSG:3857
    tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)

    shapes_3857 = []  # List of polygon geometries in EPSG:3857, clipped to tile bounds

    # Helper to transform a coordinate list
    def _reproj_coords(coords):
        for (x, y) in coords:
            yield tf.transform(x, y)

    # Clip rectangle in EPSG:3857
    clip_rect = tile_bounds_merc

    for g in candidates:
        # Fast bbox reject in lon/lat space
        if not g.bounds or not bbox_intersects(g.bounds, tile_bbox_lonlat):
            continue

        # Reproject geometry to EPSG:3857
        if isinstance(g, Polygon):
            exterior = list(_reproj_coords(g.exterior.coords))
            interiors = [list(_reproj_coords(r.coords)) for r in g.interiors]
            gg = Polygon(exterior, interiors)
        else:
            polys = []
            for p in g.geoms:
                exterior = list(_reproj_coords(p.exterior.coords))
                interiors = [list(_reproj_coords(r.coords)) for r in p.interiors]
                polys.append(Polygon(exterior, interiors))
            gg = MultiPolygon(polys)

        # Clip to tile bounds to minimize rasterization work
        gg2 = clip_by_rect(gg, clip_rect[0], clip_rect[1], clip_rect[2], clip_rect[3])
        if gg2.is_empty:
            continue

        shapes_3857.append(gg2)

    if not shapes_3857:
        return np.zeros(out_shape, dtype=bool)

    # Build a transform for the tile pixel grid
    transform = from_bounds(*tile_bounds_merc, width=out_shape[1], height=out_shape[0])

    # Rasterize polygons as 1 on land, 0 elsewhere
    mask = rasterize(
        [(geom, 1) for geom in shapes_3857],
        out_shape=out_shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
        all_touched=False,  # conservative land; set True if you prefer more masking
    )

    return mask.astype(bool)


# ============================================================
# 7) Contour extraction from raster arrays
# ============================================================

def compute_contour_levels(data: np.ndarray, interval_m: float, nodata_mask: np.ndarray) -> List[float]:
    """
    Decide contour levels based on min/max of valid data inside the tile window.
    """
    valid = data[~nodata_mask]  # Keep only non-nodata
    if valid.size == 0:
        return []

    vmin = float(np.nanmin(valid))
    vmax = float(np.nanmax(valid))

    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        return []

    # Align the level range to the interval
    lo = math.floor(vmin / interval_m) * interval_m
    hi = math.ceil(vmax / interval_m) * interval_m

    n = int(round((hi - lo) / interval_m)) + 1  # Count of levels
    if n <= 1 or n > 20000:
        return []

    return [lo + i * interval_m for i in range(n)]


def contours_from_array(
    grid: np.ndarray,
    extent_merc: Tuple[float, float, float, float],
    levels: List[float],
    nodata_mask: np.ndarray,
) -> List[Tuple[float, LineString]]:
    """
    Extract contour lines from a 2D grid using matplotlib.contour.

    Returns: list of (level_value, LineString in EPSG:3857 meters)
    """
    if not levels:
        return []

    # Make a copy and set nodata to NaN (matplotlib treats NaNs as masked)
    grid2 = grid.copy()
    grid2[nodata_mask] = np.nan

    minx, miny, maxx, maxy = extent_merc
    h, w = grid2.shape

    # Create coordinate vectors spanning the extent
    xs = np.linspace(minx, maxx, w)
    ys = np.linspace(miny, maxy, h)

    # Headless figure
    fig = plt.figure(figsize=(1, 1), dpi=72)
    ax = fig.add_subplot(111)
    ax.set_axis_off()

    # Compute contours; handle errors (e.g., all NaNs)
    try:
        cs = ax.contour(xs, ys, grid2, levels=levels, linewidths=1)
    except Exception:
        plt.close(fig)
        return []

    out: List[Tuple[float, LineString]] = []

    # Convert matplotlib paths to LineStrings
    for lev, coll in zip(cs.levels, cs.collections):
        for path in coll.get_paths():
            v = path.vertices
            if v.shape[0] < 2:
                continue
            ls = LineString(v)
            if ls.length > 0:
                out.append((float(lev), ls))

    plt.close(fig)
    return out


# ============================================================
# 8) Obfuscation utilities
# ============================================================

def add_jitter_linestring(ls: LineString, sigma_m: float, rng: random.Random) -> LineString:
    """
    Add Gaussian jitter to each vertex.
    """
    if sigma_m <= 0:
        return ls

    coords = list(ls.coords)
    jittered = [(x + rng.gauss(0.0, sigma_m), y + rng.gauss(0.0, sigma_m)) for (x, y) in coords]

    # Avoid degenerate line where first==last
    if len(jittered) >= 2 and jittered[0] == jittered[-1]:
        jittered[-1] = (jittered[-1][0] + 1e-6, jittered[-1][1] + 1e-6)

    return LineString(jittered)


def simplify_linestring(ls: LineString, tol_m: float) -> LineString:
    """
    Simplify geometry by Douglas-Peucker tolerance in meters.
    """
    if tol_m <= 0:
        return ls
    return ls.simplify(tol_m, preserve_topology=False)


# ============================================================
# 9) MVT encoding helpers
# ============================================================

def merc_to_tile_coords(
    x: float,
    y: float,
    bounds: Tuple[float, float, float, float],
    extent: int,
) -> Tuple[int, int]:
    """
    Convert EPSG:3857 coordinates into MVT tile coordinates [0..extent].
    """
    minx, miny, maxx, maxy = bounds

    # Normalize x into [0..1]
    tx = (x - minx) / (maxx - minx)

    # Normalize y into [0..1] but note tile Y increases downward
    ty = (maxy - y) / (maxy - miny)

    # Convert to integer grid
    px = int(round(tx * extent))
    py = int(round(ty * extent))

    # Clamp
    px = max(0, min(extent, px))
    py = max(0, min(extent, py))
    return px, py


def encode_mvt(layers: Dict[str, Dict[str, Any]]) -> bytes:
    """
    Encode multiple layers into a single MVT tile blob.

    layers format:
      {
        "layer_name": {
           "features": [ {geometry:{...}, properties:{...}}, ... ],
           "extent": 4096,
           "version": 2
        },
        ...
      }
    """
    payload = {"layers": []}
    for lname, layer in layers.items():
        payload["layers"].append({
            "name": lname,
            "features": layer["features"],
            "extent": layer["extent"],
            "version": layer.get("version", 2),
        })
    return mapbox_vector_tile.encode(payload)


# ============================================================
# 10) Bathymetry per-tile payload modes
# ============================================================

def _pack_grid_i16(
    data: np.ndarray,
    valid_mask: np.ndarray,
    quantization_m: float,
) -> Tuple[bytes, Dict[str, Any]]:
    """
    Pack a grid into quantized int16 array.

    - Quantization: stored_int16 = round(value / q)
    - Missing values: -32768 sentinel

    Returns:
      (compressed_bytes, metadata_properties)
    """
    q = float(quantization_m) if quantization_m and quantization_m > 0 else 1.0
    nodata_i16 = np.int16(-32768)

    out = np.full(data.shape, nodata_i16, dtype=np.int16)
    vals = np.round(data[valid_mask] / q).astype(np.int32)
    vals = np.clip(vals, -32767, 32767).astype(np.int16)
    out[valid_mask] = vals

    raw = out.tobytes(order="C")
    comp = zlib.compress(raw, level=6)

    meta = {
        "dtype": "i16",
        "q_m": q,
        "nodata": int(nodata_i16),
    }
    return comp, meta


def _pack_grid_f32(
    data: np.ndarray,
    valid_mask: np.ndarray,
) -> Tuple[bytes, Dict[str, Any]]:
    """
    Pack a grid into float32 array, using NaN for missing.

    Returns:
      (compressed_bytes, metadata_properties)
    """
    out = data.astype(np.float32, copy=True)
    out[~valid_mask] = np.nan
    raw = out.tobytes(order="C")
    comp = zlib.compress(raw, level=6)

    meta = {
        "dtype": "f32",
        "nodata": "nan",
    }
    return comp, meta


def _png_from_grid(
    data: np.ndarray,
    valid_mask: np.ndarray,
    vmin: float,
    vmax: float,
) -> bytes:
    """
    Convert grid to an 8-bit grayscale PNG in-memory.

    This is *not* a georeferenced PNG; it’s a compact per-tile payload.
    Client decodes and remaps values with provided (vmin,vmax).

    Requires Pillow at runtime.
    """
    # Import lazily so users who don't need PNG mode can avoid Pillow dependency issues
    from PIL import Image  # type: ignore

    # Avoid division by zero
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        vmin, vmax = -1.0, 1.0

    # Normalize to 0..255
    norm = (data - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)

    # Put nodata as 0 (client uses mask flag anyway)
    img = (norm * 255.0).astype(np.uint8)
    img[~valid_mask] = 0

    # Encode PNG
    im = Image.fromarray(img, mode="L")
    buf = io.BytesIO()  # noqa: F821 (io is used only if PNG mode is enabled)
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def make_bathy_payload_feature(
    src: rasterio.DatasetReader,
    tb_merc: Tuple[float, float, float, float],
    mode: str,
    grid_size: int,
    quantization_m: float,
    extent: int,
) -> Optional[Dict[str, Any]]:
    """
    Create a single MVT feature that stores per-tile bathymetry payload.

    mode options:
      - "none"         -> returns None (no bathy storage)
      - "mvt_grid_i16" -> int16 grid zlib+base64 properties
      - "mvt_grid_f32" -> float32 grid zlib+base64 properties
      - "mvt_png"      -> grayscale png zlib+base64 properties (needs Pillow)

    Returns:
      A valid MVT feature dict (geometry + properties), or None.
    """
    mode = (mode or "none").lower()
    if mode == "none":
        return None

    # Read a small grid window from the intermediate raster for this tile
    out_shape = (grid_size, grid_size)
    window = rasterio.windows.from_bounds(*tb_merc, transform=src.transform)

    data = src.read(
        1,
        window=window,
        out_shape=out_shape,
        masked=False,
        resampling=Resampling.bilinear,
    )

    nodata = src.nodata

    # Build a valid_mask (True where value is meaningful)
    if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
        valid_mask = np.isfinite(data)
    else:
        valid_mask = np.isfinite(data) & (data != nodata)

    # Prepare geometry (a single point at tile center, so the feature is legal)
    cx = (tb_merc[0] + tb_merc[2]) * 0.5
    cy = (tb_merc[1] + tb_merc[3]) * 0.5
    px, py = merc_to_tile_coords(cx, cy, tb_merc, extent)

    props: Dict[str, Any] = {
        "grid_w": int(grid_size),
        "grid_h": int(grid_size),
        "mode": mode,
    }

    # Encode payload based on mode
    if mode == "mvt_grid_i16":
        comp, meta = _pack_grid_i16(data, valid_mask, quantization_m)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")

    elif mode == "mvt_grid_f32":
        comp, meta = _pack_grid_f32(data, valid_mask)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")

    elif mode == "mvt_png":
        # PNG mode encodes values as grayscale; we also store value range for client decoding
        # Compute min/max over valid pixels
        if valid_mask.any():
            vmin = float(np.nanmin(data[valid_mask]))
            vmax = float(np.nanmax(data[valid_mask]))
        else:
            vmin, vmax = -1.0, 1.0

        # Encode PNG bytes (requires Pillow)
        import io  # Local import so it exists only when needed
        from PIL import Image  # type: ignore

        if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
            vmin, vmax = -1.0, 1.0

        norm = (data - vmin) / (vmax - vmin)
        norm = np.clip(norm, 0.0, 1.0)

        img = (norm * 255.0).astype(np.uint8)
        img[~valid_mask] = 0  # Set nodata pixels to 0

        im = Image.fromarray(img, mode="L")
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=True)
        png_bytes = buf.getvalue()

        # Compress PNG bytes (often smaller after zlib) then base64
        comp = zlib.compress(png_bytes, level=6)
        props["dtype"] = "png8"
        props["vmin"] = vmin
        props["vmax"] = vmax
        props["nodata_is_zero"] = True
        props["png_b64"] = base64.b64encode(comp).decode("ascii")

    else:
        raise ValueError(f"Unknown bathy storage mode: {mode}")

    return {
        "geometry": {"type": "Point", "coordinates": [px, py]},
        "properties": props,
    }


# ============================================================
# 11) Mosaic + per-zoom intermediate cache
# ============================================================

def build_mosaic_sources(tile_tifs: List[str]) -> List[rasterio.DatasetReader]:
    """
    Open all tile GeoTIFFs for merging.
    """
    return [rasterio.open(p) for p in tile_tifs]


def close_mosaic_sources(srcs: List[rasterio.DatasetReader]) -> None:
    """
    Close rasterio datasets.
    """
    for s in srcs:
        try:
            s.close()
        except Exception:
            pass


def ensure_zoom_intermediate(
    cfg: "Config",
    tile_tifs: List[str],
    z: int,
    bbox_merc: Tuple[float, float, float, float],
) -> str:
    """
    Build or reuse a cached intermediate raster for zoom z:
      - mosaics the source tiles
      - warps to EPSG:3857
      - clips to bbox
      - resamples to zoom-appropriate resolution

    This drastically speeds up per-tile reads.
    """
    os.makedirs(cfg.cache_dir, exist_ok=True)
    z_dir = os.path.join(cfg.cache_dir, f"z{z}")
    os.makedirs(z_dir, exist_ok=True)

    out_tif = os.path.join(z_dir, f"dtm_{cfg.dtm_year}_bbox_z{z}_3857.tif")

    # If caching enabled and file exists, reuse it
    if cfg.cache_intermediate and os.path.exists(out_tif):
        return out_tif

    # Mosaic in source CRS
    srcs = build_mosaic_sources(tile_tifs)
    try:
        mosaic, mosaic_transform = rio_merge(srcs)  # (bands, H, W)
        data = mosaic[0]  # Take band 1
        src_crs = srcs[0].crs
        nodata = srcs[0].nodata
        if nodata is None:
            nodata = np.nan
    finally:
        close_mosaic_sources(srcs)

    # Build destination grid in EPSG:3857
    minx, miny, maxx, maxy = bbox_merc
    res = zoom_target_resolution_m(z, cfg.tile_pixels, tile_size=256)

    width = max(1, int(math.ceil((maxx - minx) / res)))
    height = max(1, int(math.ceil((maxy - miny) / res)))

    dst_transform = from_bounds(minx, miny, maxx, maxy, width=width, height=height)

    # Allocate destination array filled with nodata
    dst = np.full((height, width), nodata, dtype=data.dtype)

    # Reproject mosaic -> EPSG:3857
    rasterio.warp.reproject(
        source=data,
        destination=dst,
        src_transform=mosaic_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs="EPSG:3857",
        resampling=Resampling.bilinear,
        src_nodata=nodata,
        dst_nodata=nodata,
    )

    # Write cached GeoTIFF (tiled + DEFLATE for reasonable disk use and performance)
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": dst.dtype,
        "crs": "EPSG:3857",
        "transform": dst_transform,
        "nodata": nodata,
        "compress": "DEFLATE",
        "predictor": 2,
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(out_tif, "w", **profile) as w:
        w.write(dst, 1)

    return out_tif


# ============================================================
# 12) Configuration model
# ============================================================

@dataclass
class Config:
    # EMODnet discovery
    dtm_year: int
    csw_url: str

    # Paths
    output_mbtiles: str
    cache_dir: str

    # Area & zoom
    bbox_lonlat: Tuple[float, float, float, float]
    zoom_min: int
    zoom_max: int

    # Vector tiling
    tile_pixels: int
    extent: int
    contour_layer_name: str
    contour_intervals_m: Dict[str, float]

    # Intermediate caching
    cache_intermediate: bool

    # Coastline mask
    coast_enabled: bool
    coast_source: str
    osm_land_polygons_url: str
    coast_max_features: Optional[int]

    # Per-tile bathy payload (client-side processing)
    bathy_enabled: bool
    bathy_mode: str            # none | mvt_grid_i16 | mvt_grid_f32 | mvt_png
    bathy_layer_name: str
    bathy_grid_size: int
    bathy_quantization_m: float

    # Obfuscation
    ob_enabled: bool
    ob_seed: int
    ob_jitter_m: Dict[str, float]
    ob_simplify_m: Dict[str, float]
    ob_drop_p: Dict[str, float]
    ob_level_noise_m: Dict[str, float]

    # MBTiles metadata
    metadata: Dict[str, str]


def load_config(path: str) -> Config:
    """
    Load JSON config and normalize to Config dataclass.
    """
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    def req(k: str):
        if k not in raw:
            raise ValueError(f"Missing required config key: {k}")
        return raw[k]

    ob = raw.get("obfuscation", {})
    coast = raw.get("coastline_mask", {})
    bathy = raw.get("bathy_data", {})  # NEW: unified bathy payload options

    return Config(
        dtm_year=int(raw.get("dtm_year", 2024)),
        csw_url=str(raw.get("csw_url", "https://emodnet.ec.europa.eu/geonetwork/emodnet/eng/csw")),

        output_mbtiles=req("output_mbtiles"),
        cache_dir=str(raw.get("cache_dir", "./cache")),

        bbox_lonlat=tuple(req("bbox_lonlat")),
        zoom_min=int(req("zoom_min")),
        zoom_max=int(req("zoom_max")),

        tile_pixels=int(raw.get("tile_pixels", 512)),
        extent=int(raw.get("extent", 4096)),
        contour_layer_name=str(raw.get("layer_name", "bathy_contours")),
        contour_intervals_m=dict(req("contour_intervals_m")),

        cache_intermediate=bool(raw.get("cache_intermediate", True)),

        coast_enabled=bool(coast.get("enabled", True)),
        coast_source=str(coast.get("source", "osm_land_polygons")),
        osm_land_polygons_url=str(coast.get(
            "osm_land_polygons_url",
            "https://osmdata.openstreetmap.de/data/land-polygons-split-4326.zip"
        )),
        coast_max_features=(int(coast["max_features"]) if "max_features" in coast and coast["max_features"] is not None else None),

        bathy_enabled=bool(bathy.get("enabled", False)),
        bathy_mode=str(bathy.get("mode", "none")),
        bathy_layer_name=str(bathy.get("layer_name", "bathy_data")),
        bathy_grid_size=int(bathy.get("grid_size", 32)),
        bathy_quantization_m=float(bathy.get("quantization_m", 0.1)),

        ob_enabled=bool(ob.get("enabled", True)),
        ob_seed=int(ob.get("seed", 12345)),
        ob_jitter_m=dict(ob.get("jitter_m_by_zoom", {})),
        ob_simplify_m=dict(ob.get("simplify_m_by_zoom", {})),
        ob_drop_p=dict(ob.get("drop_probability_by_zoom", {})),
        ob_level_noise_m=dict(ob.get("level_noise_m_by_zoom", {})),

        metadata=dict(raw.get("mbtiles_metadata", {})),
    )


# ============================================================
# 13) Tile iteration
# ============================================================

def tile_iter_for_bbox(bbox_lonlat: Tuple[float, float, float, float], z: int) -> Iterable[mercantile.Tile]:
    """
    Yield mercantile.Tile objects for all tiles intersecting bbox at zoom z.
    """
    minlon, minlat, maxlon, maxlat = bbox_lonlat
    return mercantile.tiles(minlon, minlat, maxlon, maxlat, [z])


# ============================================================
# 14) Main pipeline
# ============================================================

def run(cfg: Config) -> None:
    """
    Orchestrate the entire pipeline using config.
    """
    rng = random.Random(cfg.ob_seed)  # Deterministic randomness for reproducible obfuscation

    bbox_merc = lonlat_bbox_to_merc(cfg.bbox_lonlat)  # EPSG:3857 bbox for fast intersection checks
    bbox_geom_merc = box(*bbox_merc)  # Shapely bbox polygon

    # ------------------------
    # A) Discover EMODnet tiles
    # ------------------------
    print("Discovering EMODnet DTM tiles via CSW...")
    tile_urls = csw_search_tiles(cfg.csw_url, cfg.bbox_lonlat, cfg.dtm_year)
    if not tile_urls:
        raise RuntimeError("No EMODnet DTM tile URLs found for your bbox/year. Check bbox, dtm_year, csw_url.")

    print(f"Found {len(tile_urls)} tile downloads.")

    # ------------------------
    # B) Download tiles (cached)
    # ------------------------
    tiles_dir = os.path.join(cfg.cache_dir, "tiles", f"dtm_{cfg.dtm_year}")
    os.makedirs(tiles_dir, exist_ok=True)

    tile_tifs: List[str] = []
    for u in tile_urls:
        print(f"Downloading/extracting: {u}")
        tif_path = download_and_extract_zip(u, tiles_dir)
        tile_tifs.append(tif_path)

    # ------------------------
    # C) Coastline/land mask setup
    # ------------------------
    land_tree: Optional[STRtree] = None
    land_geoms: List[Polygon | MultiPolygon] = []

    if cfg.coast_enabled:
        if cfg.coast_source != "osm_land_polygons":
            raise ValueError("Currently supported coastline_mask.source: osm_land_polygons")

        print("Preparing coastline/land mask from OSM land polygons (cached)...")
        land_root = ensure_osm_land_polygons(cfg.cache_dir, cfg.osm_land_polygons_url)
        land_tree, land_geoms = load_land_polygons_strtree(
            land_root,
            cfg.bbox_lonlat,
            max_features=cfg.coast_max_features,
        )
        print(f"Loaded land polygons: {len(land_geoms)} (bbox-filtered)")

    # ------------------------
    # D) Open MBTiles + metadata
    # ------------------------
    conn = mbtiles_open(cfg.output_mbtiles)

    md = {
        "name": cfg.metadata.get("name", "Bathymetry contours (EMODnet-derived)"),
        "format": "pbf",
        "minzoom": str(cfg.zoom_min),
        "maxzoom": str(cfg.zoom_max),
        "bounds": ",".join(map(str, cfg.bbox_lonlat)),
        "attribution": cfg.metadata.get("attribution", "Derived from EMODnet Bathymetry DTM (CC BY 4.0)."),
    }
    md.update(cfg.metadata)
    mbtiles_put_metadata(conn, md)

    # ------------------------
    # E) Process zoom levels
    # ------------------------
    total_tiles = 0
    written_tiles = 0
    t0 = time.time()

    for z in range(cfg.zoom_min, cfg.zoom_max + 1):
        interval = _value_for_zoom(cfg.contour_intervals_m, z, None)
        if interval is None or interval <= 0:
            print(f"[z={z}] skip (no contour interval configured)")
            continue

        # Create/reuse intermediate cached raster for this zoom
        inter_tif = ensure_zoom_intermediate(cfg, tile_tifs, z, bbox_merc)

        # Read obfuscation parameters (used only if cfg.ob_enabled == True)
        sigma_jitter = float(_value_for_zoom(cfg.ob_jitter_m, z, 0.0) or 0.0)
        tol_simplify = float(_value_for_zoom(cfg.ob_simplify_m, z, 0.0) or 0.0)
        drop_p = float(_value_for_zoom(cfg.ob_drop_p, z, 0.0) or 0.0)
        level_noise = float(_value_for_zoom(cfg.ob_level_noise_m, z, 0.0) or 0.0)

        print(
            f"[z={z}] interval={interval}m intermediate={os.path.basename(inter_tif)} "
            f"obfuscation={'ON' if cfg.ob_enabled else 'OFF'} "
            f"land_mask={'ON' if cfg.coast_enabled else 'OFF'} "
            f"bathy_data={'ON' if cfg.bathy_enabled else 'OFF'}({cfg.bathy_mode})"
        )

        tiles = list(tile_iter_for_bbox(cfg.bbox_lonlat, z))
        total_tiles += len(tiles)

        # Open intermediate raster once per zoom for efficiency
        with rasterio.open(inter_tif) as src:
            nodata = src.nodata

            for t in tiles:
                x, y = t.x, t.y
                tb = tile_bounds_merc(z, x, y)  # tile bounds in EPSG:3857
                tile_geom = box(*tb)

                # Skip tiles fully outside bbox
                if not tile_geom.intersects(bbox_geom_merc):
                    continue

                # Read bathy window for contouring at tile_pixels resolution
                out_shape = (cfg.tile_pixels, cfg.tile_pixels)
                window = rasterio.windows.from_bounds(*tb, transform=src.transform)

                try:
                    data = src.read(
                        1,
                        window=window,
                        out_shape=out_shape,
                        masked=False,
                        resampling=Resampling.bilinear,
                    )
                except Exception:
                    # Window might go outside raster, or read fails for other reasons
                    continue

                # Build nodata mask
                if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                    nodata_mask = ~np.isfinite(data)
                else:
                    nodata_mask = (data == nodata) | ~np.isfinite(data)

                # Apply land mask: treat land pixels as nodata, so contours will not be generated there
                if cfg.coast_enabled and land_tree is not None:
                    land_mask = rasterize_land_mask(land_tree, land_geoms, tb, out_shape)
                    if land_mask.any():
                        nodata_mask = nodata_mask | land_mask

                # ------------------------
                # Contour extraction
                # ------------------------
                levels = compute_contour_levels(data, float(interval), nodata_mask)
                raw_contours = contours_from_array(data, tb, levels, nodata_mask) if levels else []

                contour_features: List[Dict[str, Any]] = []

                # Convert contours to MVT features
                for lev, ls in raw_contours:
                    # Clip strictly to tile bounds
                    ls2 = clip_by_rect(ls, tb[0], tb[1], tb[2], tb[3])
                    if ls2.is_empty:
                        continue

                    # MultiLineString support
                    lines = [ls2] if isinstance(ls2, LineString) else list(ls2.geoms)

                    for line in lines:
                        if line.length <= 0:
                            continue

                        # Optional obfuscation
                        if cfg.ob_enabled:
                            # Random feature drop
                            if drop_p > 0 and rng.random() < drop_p:
                                continue

                            # Jitter vertices
                            line = add_jitter_linestring(line, sigma_jitter, rng)

                            # Simplify geometry
                            line = simplify_linestring(line, tol_simplify)

                            # Discard degenerate lines
                            if line.is_empty or len(line.coords) < 2:
                                continue

                        # Optional level (depth) noise
                        lev2 = float(lev)
                        if cfg.ob_enabled and level_noise > 0:
                            lev2 += rng.gauss(0.0, level_noise)

                        # Convert geometry to MVT LineString coordinates
                        geom = {
                            "type": "LineString",
                            "coordinates": [merc_to_tile_coords(xx, yy, tb, cfg.extent) for (xx, yy) in line.coords],
                        }

                        contour_features.append({
                            "geometry": geom,
                            "properties": {
                                "depth_m": round(lev2, 3),
                                "interval_m": float(interval),
                                "z": int(z),
                            },
                        })

                # ------------------------
                # Per-tile bathy payload feature (optional)
                # ------------------------
                bathy_feature: Optional[Dict[str, Any]] = None
                if cfg.bathy_enabled:
                    try:
                        bathy_feature = make_bathy_payload_feature(
                            src=src,
                            tb_merc=tb,
                            mode=cfg.bathy_mode,
                            grid_size=cfg.bathy_grid_size,
                            quantization_m=cfg.bathy_quantization_m,
                            extent=cfg.extent,
                        )
                    except Exception as e:
                        # If bathy payload fails, continue with contours only
                        bathy_feature = None

                # If tile has no content at all, skip writing it
                if not contour_features and bathy_feature is None:
                    continue

                # Build layers for this tile
                layers: Dict[str, Dict[str, Any]] = {}

                if contour_features:
                    layers[cfg.contour_layer_name] = {
                        "features": contour_features,
                        "extent": cfg.extent,
                        "version": 2,
                    }

                if bathy_feature is not None:
                    layers[cfg.bathy_layer_name] = {
                        "features": [bathy_feature],
                        "extent": cfg.extent,
                        "version": 2,
                    }

                # Encode MVT tile blob
                tile_data = encode_mvt(layers)

                # Write to MBTiles
                mbtiles_put_tile(conn, z, x, tms_y(z, y), tile_data)
                written_tiles += 1

        # Commit per zoom (prevents huge transactions)
        mbtiles_commit(conn)

    # ------------------------
    # F) Wrap up
    # ------------------------
    dt = time.time() - t0
    print(f"Done. total_tiles={total_tiles}, written_tiles={written_tiles}, elapsed={dt:.1f}s")

    conn.commit()
    conn.close()


# ============================================================
# 15) CLI entrypoint
# ============================================================

def main(argv: List[str]) -> int:
    """
    Parse CLI arguments and run pipeline.
    """
    if len(argv) != 2:
        print("Usage: python obfuscate_bathy_tiles.py <config.json>", file=sys.stderr)
        return 2

    cfg = load_config(argv[1])
    run(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
