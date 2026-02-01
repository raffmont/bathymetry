#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
obfuscate_bathy_tiles.py

Generate bathymetric contour vector tiles (MVT inside MBTiles) from EMODnet DTM tiles,
with optional obfuscation and land masking.

New in this version
-------------------
1) Bathymetry payload output can be:
   - embedded in the vector tile (MVT properties), OR
   - written to a separate raster MBTiles as PNG tiles, OR
   - disabled entirely.

2) Land masking can be sourced from:
   - OSM land polygons dataset (downloaded & cached), OR
   - an EMODnet coastline product you provide as a direct vector download URL
     (ZIP shapefile, GeoPackage, GeoJSON, etc.).

3) Land mask is cached per zoom as a GeoTIFF aligned to the per-zoom intermediate raster,
   so per-tile masking is a fast window read instead of per-tile polygon rasterization.

Usage:
  python obfuscate_bathy_tiles.py config.json

Dependencies:
  See requirements.txt from prior step.
"""

from __future__ import annotations

import base64
import io
import json
import math
import os
import random
import re
import shutil
import sqlite3
import sys
import time
import zipfile
import zlib
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import mercantile
import numpy as np
import requests
import rasterio
import rasterio.warp
from pyproj import CRS, Transformer
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.merge import merge as rio_merge
from rasterio.transform import from_bounds

import fiona
from shapely.geometry import LineString, Polygon, MultiPolygon, box, shape, mapping
from shapely.ops import clip_by_rect
from shapely.strtree import STRtree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mapbox_vector_tile

# Optional dependency for PNG raster MBTiles / embedded png payload
try:
    from PIL import Image
    PIL_OK = True
except Exception:
    PIL_OK = False


# =========================
# Zoom-range config helpers
# =========================

def _parse_zoom_key(key: str) -> Tuple[int, int]:
    key = key.strip()
    if "-" in key:
        a, b = key.split("-", 1)
        return int(a), int(b)
    z = int(key)
    return z, z

def _value_for_zoom(ranged: Dict[str, Any], z: int, default: Any = None) -> Any:
    exact = str(z)
    if exact in ranged:
        return ranged[exact]
    for k, v in ranged.items():
        z0, z1 = _parse_zoom_key(k)
        if z0 <= z <= z1:
            return v
    return default


# ==============
# MBTiles schema
# ==============

MBTILES_SCHEMA = """
CREATE TABLE IF NOT EXISTS metadata (name TEXT, value TEXT);
CREATE TABLE IF NOT EXISTS tiles (zoom_level INTEGER, tile_column INTEGER, tile_row INTEGER, tile_data BLOB);
CREATE UNIQUE INDEX IF NOT EXISTS tile_index ON tiles (zoom_level, tile_column, tile_row);
"""

def mbtiles_open(path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.executescript(MBTILES_SCHEMA)
    return conn

def mbtiles_put_metadata(conn: sqlite3.Connection, md: Dict[str, str]) -> None:
    conn.execute("DELETE FROM metadata;")
    for k, v in md.items():
        conn.execute("INSERT INTO metadata(name, value) VALUES(?, ?);", (k, v))
    conn.commit()

def mbtiles_put_tile(conn: sqlite3.Connection, z: int, x: int, y_tms: int, data: bytes) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO tiles(zoom_level, tile_column, tile_row, tile_data) VALUES(?,?,?,?);",
        (z, x, y_tms, data),
    )

def mbtiles_commit(conn: sqlite3.Connection) -> None:
    conn.commit()

def tms_y(z: int, y_xyz: int) -> int:
    return (2**z - 1) - y_xyz


# =====================
# WebMercator helpers
# =====================

WEBMERC_MAX = 20037508.342789244

def lonlat_bbox_to_merc(b: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    minlon, minlat, maxlon, maxlat = b
    tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    minx, miny = tf.transform(minlon, minlat)
    maxx, maxy = tf.transform(maxlon, maxlat)
    return (minx, miny, maxx, maxy)

def merc_bbox_to_lonlat(b: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    minx, miny, maxx, maxy = b
    tf = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    minlon, minlat = tf.transform(minx, miny)
    maxlon, maxlat = tf.transform(maxx, maxy)
    return (min(minlon, maxlon), min(minlat, maxlat), max(minlon, maxlon), max(minlat, maxlat))

def tile_bounds_merc(z: int, x: int, y: int) -> Tuple[float, float, float, float]:
    b = mercantile.bounds(x, y, z)  # lon/lat
    return lonlat_bbox_to_merc((b.west, b.south, b.east, b.north))

def meters_per_pixel_equator(z: int, tile_size: int = 256) -> float:
    world = 2 * WEBMERC_MAX
    return world / (tile_size * 2**z)

def zoom_target_resolution_m(z: int, tile_pixels: int, tile_size: int = 256) -> float:
    return meters_per_pixel_equator(z, tile_size=tile_size) * (tile_size / tile_pixels)


# ======================
# EMODnet CSW discovery
# ======================

CSW_GETRECORDS_TEMPLATE = """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<csw:GetRecords
  xmlns:csw=\"http://www.opengis.net/cat/csw/2.0.2\"
  xmlns:ogc=\"http://www.opengis.net/ogc\"
  service=\"CSW\"
  version=\"2.0.2\"
  resultType=\"results\"
  startPosition=\"{start}\"
  maxRecords=\"{max_records}\"
  outputSchema=\"http://www.isotc211.org/2005/gmd\"
  outputFormat=\"application/xml\">
  <csw:Query typeNames=\"csw:Record\">
    <csw:ElementSetName>full</csw:ElementSetName>
    <csw:Constraint version=\"1.1.0\">
      <ogc:Filter>
        <ogc:And>
          <ogc:PropertyIsLike wildCard=\"*\" singleChar=\"?\" escapeChar=\"\\\">
            <ogc:PropertyName>AnyText</ogc:PropertyName>
            <ogc:Literal>*EMODnet Digital Bathymetry (DTM {year}) - Tile*</ogc:Literal>
          </ogc:PropertyIsLike>
          <ogc:BBOX>
            <ogc:PropertyName>ows:BoundingBox</ogc:PropertyName>
            <gml:Envelope xmlns:gml=\"http://www.opengis.net/gml\">
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
    minx, miny, maxx, maxy = bbox_lonlat
    tile_urls: List[str] = []
    seen = set()

    url_re = re.compile(r"https://downloads\.emodnet-bathymetry\.eu/[^\s\"<]+?_{}\.tif\.zip".format(year))

    start = 1
    max_records = 50
    headers = {"Content-Type": "application/xml"}

    for _ in range(max_pages):
        body = CSW_GETRECORDS_TEMPLATE.format(
            start=start,
            max_records=max_records,
            year=year,
            minx=minx,
            miny=miny,
            maxx=maxx,
            maxy=maxy,
        )
        r = requests.post(csw_url, data=body.encode("utf-8"), headers=headers, timeout=timeout_s)
        r.raise_for_status()
        xml = r.text

        found = url_re.findall(xml)
        for u in found:
            if u not in seen:
                seen.add(u)
                tile_urls.append(u)

        if "nextRecord=\"0\"" in xml or "numberOfRecordsReturned=\"0\"" in xml:
            break
        if len(found) == 0 and start > 1:
            break
        start += max_records

    return tile_urls


# ======================
# Download & extraction
# ======================

def download_file(url: str, out_path: str, timeout_s: int = 180) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    if os.path.exists(out_path):
        return
    with requests.get(url, stream=True, timeout=timeout_s) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

def download_and_extract_zip(url: str, out_dir: str, timeout_s: int = 180) -> str:
    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.basename(url)
    local_zip = os.path.join(out_dir, fname)

    download_file(url, local_zip, timeout_s=timeout_s)

    with zipfile.ZipFile(local_zip, "r") as z:
        tif_members = [m for m in z.namelist() if m.lower().endswith((".tif", ".tiff"))]
        if not tif_members:
            raise RuntimeError(f"No GeoTIFF found inside {local_zip}")
        member = tif_members[0]
        out_tif = os.path.join(out_dir, os.path.basename(member))
        if not os.path.exists(out_tif):
            z.extract(member, out_dir)
            extracted = os.path.join(out_dir, member)
            if extracted != out_tif:
                os.makedirs(os.path.dirname(out_tif), exist_ok=True)
                shutil.move(extracted, out_tif)
        return out_tif

def download_and_prepare_vector_dataset(url: str, cache_dir: str) -> Tuple[str, Optional[str]]:
    """
    Download a vector dataset to cache and return (dataset_path, layer_name_or_None).
    Supports:
      - ZIP (shapefile bundle or other Fiona-readable)
      - direct .gpkg / .geojson / .json / .shp (and friends)
    If ZIP: extracted into cache folder; returns extracted root path; layer None (Fiona will list layers).
    """
    os.makedirs(cache_dir, exist_ok=True)
    fname = os.path.basename(url.split("?", 1)[0])
    local_path = os.path.join(cache_dir, fname)

    download_file(url, local_path, timeout_s=600)

    if fname.lower().endswith(".zip"):
        extract_dir = os.path.join(cache_dir, fname[:-4])
        marker = os.path.join(extract_dir, ".extracted")
        if not os.path.exists(marker):
            os.makedirs(extract_dir, exist_ok=True)
            with zipfile.ZipFile(local_path, "r") as z:
                z.extractall(extract_dir)
            with open(marker, "w", encoding="utf-8") as f:
                f.write("ok\n")
        return extract_dir, None

    return local_path, None

def find_vector_files(root_dir: str) -> List[str]:
    """Find candidate vector files under extracted zip root."""
    exts = (".shp", ".gpkg", ".geojson", ".json")
    out = []
    for base, _, files in os.walk(root_dir):
        for fn in files:
            if fn.lower().endswith(exts):
                out.append(os.path.join(base, fn))
    return out


# ============================
# Coastline masking data load
# ============================

def ensure_osm_land_polygons(cache_dir: str, url: str) -> str:
    coast_dir = os.path.join(cache_dir, "coastline_mask", "osm_land_polygons")
    os.makedirs(coast_dir, exist_ok=True)
    zip_path = os.path.join(coast_dir, os.path.basename(url))
    download_file(url, zip_path, timeout_s=600)

    marker = os.path.join(coast_dir, ".extracted")
    if not os.path.exists(marker):
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(coast_dir)
        with open(marker, "w", encoding="utf-8") as f:
            f.write("ok\n")
    return coast_dir

def bbox_intersects(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> bool:
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])

def _transform_geom(geom, tf: Transformer):
    """Reproject shapely geometry by transforming all coordinates."""
    if geom.is_empty:
        return geom
    gtype = geom.geom_type
    if gtype == "Polygon":
        ext = [tf.transform(x, y) for (x, y) in geom.exterior.coords]
        holes = [[tf.transform(x, y) for (x, y) in r.coords] for r in geom.interiors]
        return Polygon(ext, holes)
    if gtype == "MultiPolygon":
        polys = []
        for p in geom.geoms:
            ext = [tf.transform(x, y) for (x, y) in p.exterior.coords]
            holes = [[tf.transform(x, y) for (x, y) in r.coords] for r in p.interiors]
            polys.append(Polygon(ext, holes))
        return MultiPolygon(polys)
    # Fallback: use mapping -> coords recursion (slower) not implemented; skip non-polys
    return geom

def load_polygons_from_fiona_source(
    path: str,
    bbox_lonlat: Tuple[float, float, float, float],
    layer: Optional[str] = None,
    max_features: Optional[int] = None,
) -> List[Polygon | MultiPolygon]:
    """
    Load polygonal geometries from a Fiona-readable dataset, filtering by bbox if possible.
    Normalizes output geometries into EPSG:4326 so we can build one STRtree.
    """
    geoms: List[Polygon | MultiPolygon] = []
    count = 0

    # If container has multiple layers, choose one
    with fiona.open(path, layer=layer) as src:
        src_crs = CRS.from_user_input(src.crs_wkt or src.crs or "EPSG:4326")
        dst_crs = CRS.from_epsg(4326)

        # Transformer from src CRS -> EPSG:4326 (lon/lat)
        tf_to_4326 = Transformer.from_crs(src_crs, dst_crs, always_xy=True)

        # Attempt bbox filter in source CRS by transforming bbox corners
        try:
            if src_crs != dst_crs:
                tf_bbox = Transformer.from_crs(dst_crs, src_crs, always_xy=True)
                minlon, minlat, maxlon, maxlat = bbox_lonlat
                x1, y1 = tf_bbox.transform(minlon, minlat)
                x2, y2 = tf_bbox.transform(maxlon, maxlat)
                bbox_src = (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
            else:
                bbox_src = bbox_lonlat
        except Exception:
            bbox_src = None

        # Iterate features
        for feat in src:
            g = feat.get("geometry")
            if not g:
                continue
            s = shape(g)
            if s.is_empty:
                continue
            if not isinstance(s, (Polygon, MultiPolygon)):
                continue

            # Fast reject by bbox in source CRS if we could compute it
            if bbox_src is not None:
                gb = s.bounds
                if gb and not bbox_intersects(gb, bbox_src):
                    continue

            # Reproject into EPSG:4326
            s4326 = _transform_geom(s, tf_to_4326)

            # Final reject against lon/lat bbox
            if not bbox_intersects(s4326.bounds, bbox_lonlat):
                continue

            geoms.append(s4326)
            count += 1
            if max_features and count >= max_features:
                break

    return geoms

def load_land_mask_source(cfg, bbox_lonlat: Tuple[float, float, float, float]) -> Tuple[STRtree, List[Polygon | MultiPolygon]]:
    """Load land/coast polygons in EPSG:4326 and build STRtree."""
    source = cfg.coast_source

    if source == "osm_land_polygons":
        root = ensure_osm_land_polygons(cfg.cache_dir, cfg.osm_land_polygons_url)
        # OSM land polygons zip contains multiple shapefiles; load only those intersecting bbox by bounds
        shps = []
        for base, _, files in os.walk(root):
            for fn in files:
                if fn.lower().endswith(".shp"):
                    shps.append(os.path.join(base, fn))

        geoms: List[Polygon | MultiPolygon] = []
        for shp in shps:
            # Quick skip by shapefile bounds (EPSG:4326 in this dataset)
            try:
                with fiona.open(shp) as src:
                    if not bbox_intersects(src.bounds, bbox_lonlat):
                        continue
            except Exception:
                continue
            geoms.extend(load_polygons_from_fiona_source(shp, bbox_lonlat, layer=None, max_features=cfg.coast_max_features))
            if cfg.coast_max_features and len(geoms) >= cfg.coast_max_features:
                geoms = geoms[: cfg.coast_max_features]
                break

        tree = STRtree(geoms) if geoms else STRtree([])
        return tree, geoms

    if source == "emodnet_product":
        # Download dataset (zip or direct) into cache
        vcache = os.path.join(cfg.cache_dir, "coastline_mask", "emodnet_product")
        ds_path, _ = download_and_prepare_vector_dataset(cfg.emodnet_vector_url, vcache)

        # If zip extracted, we need to choose a file
        if os.path.isdir(ds_path):
            candidates = find_vector_files(ds_path)
            if not candidates:
                raise RuntimeError(f"No vector files found in extracted EMODnet dataset: {ds_path}")
            # Prefer .gpkg, else .shp, else .geojson
            candidates.sort(key=lambda p: (0 if p.lower().endswith(".gpkg") else 1 if p.lower().endswith(".shp") else 2, p))
            ds_path = candidates[0]

        geoms = load_polygons_from_fiona_source(
            ds_path,
            bbox_lonlat,
            layer=cfg.emodnet_vector_layer,
            max_features=cfg.coast_max_features,
        )
        tree = STRtree(geoms) if geoms else STRtree([])
        return tree, geoms

    raise ValueError(f"Unsupported coastline_mask.source: {source}")


# =========================================
# Land mask raster cache (per zoom GeoTIFF)
# =========================================

def ensure_landmask_zoom(
    cfg,
    z: int,
    inter_tif: str,
    bbox_merc: Tuple[float, float, float, float],
    land_tree: STRtree,
    land_geoms: List[Polygon | MultiPolygon],
) -> str:
    """Create or reuse a per-zoom landmask GeoTIFF aligned to the intermediate raster."""
    z_dir = os.path.join(cfg.cache_dir, f"z{z}")
    os.makedirs(z_dir, exist_ok=True)
    out_mask = os.path.join(z_dir, f"landmask_z{z}_3857.tif")

    if cfg.landmask_cache_enabled and os.path.exists(out_mask):
        return out_mask

    # If no land geometries, write all-zero mask
    with rasterio.open(inter_tif) as src:
        profile = src.profile.copy()
        profile.update({
            "count": 1,
            "dtype": "uint8",
            "nodata": 0,
            "compress": "DEFLATE",
            "predictor": 2,
        })
        H, W = src.height, src.width
        transform = src.transform

        if not land_geoms:
            mask_arr = np.zeros((H, W), dtype=np.uint8)
            with rasterio.open(out_mask, "w", **profile) as dst:
                dst.write(mask_arr, 1)
            return out_mask

        # Query relevant polygons in lon/lat bbox (same as config bbox)
        bbox_lonlat = cfg.bbox_lonlat
        bbox_poly = box(*bbox_lonlat)
        candidates = land_tree.query(bbox_poly)

        # Reproject candidates to EPSG:3857 and clip to bbox_merc
        tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
        minx, miny, maxx, maxy = bbox_merc
        shapes_3857 = []
        for g in candidates:
            if not bbox_intersects(g.bounds, bbox_lonlat):
                continue
            gg = _transform_geom(g, tf)
            gg = clip_by_rect(gg, minx, miny, maxx, maxy)
            if gg.is_empty:
                continue
            shapes_3857.append(gg)

        # Rasterize on the intermediate raster grid
        mask_arr = rasterize(
            [(geom, 1) for geom in shapes_3857],
            out_shape=(H, W),
            transform=transform,
            fill=0,
            dtype=np.uint8,
            all_touched=bool(cfg.landmask_all_touched),
        )

        with rasterio.open(out_mask, "w", **profile) as dst:
            dst.write(mask_arr, 1)

    return out_mask


# ==================================
# Contouring & MVT encoding helpers
# ==================================

def compute_contour_levels(data: np.ndarray, interval_m: float, nodata_mask: np.ndarray) -> List[float]:
    valid = data[~nodata_mask]
    if valid.size == 0:
        return []
    vmin = float(np.nanmin(valid))
    vmax = float(np.nanmax(valid))
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        return []
    lo = math.floor(vmin / interval_m) * interval_m
    hi = math.ceil(vmax / interval_m) * interval_m
    n = int(round((hi - lo) / interval_m)) + 1
    if n <= 1 or n > 20000:
        return []
    return [lo + i * interval_m for i in range(n)]

def contours_from_array(
    grid: np.ndarray,
    extent_merc: Tuple[float, float, float, float],
    levels: List[float],
    nodata_mask: np.ndarray,
) -> List[Tuple[float, LineString]]:
    if not levels:
        return []
    grid2 = grid.copy()
    grid2[nodata_mask] = np.nan

    minx, miny, maxx, maxy = extent_merc
    h, w = grid2.shape
    xs = np.linspace(minx, maxx, w)
    ys = np.linspace(miny, maxy, h)

    fig = plt.figure(figsize=(1, 1), dpi=72)
    ax = fig.add_subplot(111)
    ax.set_axis_off()

    try:
        cs = ax.contour(xs, ys, grid2, levels=levels, linewidths=1)
    except Exception:
        plt.close(fig)
        return []

    out: List[Tuple[float, LineString]] = []
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

def add_jitter_linestring(ls: LineString, sigma_m: float, rng: random.Random) -> LineString:
    if sigma_m <= 0:
        return ls
    coords = list(ls.coords)
    jittered = [(x + rng.gauss(0.0, sigma_m), y + rng.gauss(0.0, sigma_m)) for (x, y) in coords]
    if len(jittered) >= 2 and jittered[0] == jittered[-1]:
        jittered[-1] = (jittered[-1][0] + 1e-6, jittered[-1][1] + 1e-6)
    return LineString(jittered)

def simplify_linestring(ls: LineString, tol_m: float) -> LineString:
    if tol_m <= 0:
        return ls
    return ls.simplify(tol_m, preserve_topology=False)

def merc_to_tile_coords(
    x: float, y: float,
    bounds: Tuple[float, float, float, float],
    extent: int
) -> Tuple[int, int]:
    minx, miny, maxx, maxy = bounds
    tx = (x - minx) / (maxx - minx)
    ty = (maxy - y) / (maxy - miny)
    px = int(round(tx * extent))
    py = int(round(ty * extent))
    px = max(0, min(extent, px))
    py = max(0, min(extent, py))
    return px, py

def encode_mvt(layers: Dict[str, Dict[str, Any]]) -> bytes:
    payload = {"layers": []}
    for lname, layer in layers.items():
        payload["layers"].append({
            "name": lname,
            "features": layer["features"],
            "extent": layer["extent"],
            "version": layer.get("version", 2),
        })
    return mapbox_vector_tile.encode(payload)


# ============================
# Bathy payload: embedded MVT
# ============================

def _pack_grid_i16(data: np.ndarray, valid_mask: np.ndarray, q_m: float) -> Tuple[bytes, Dict[str, Any]]:
    q = float(q_m) if q_m and q_m > 0 else 1.0
    nodata_i16 = np.int16(-32768)
    out = np.full(data.shape, nodata_i16, dtype=np.int16)
    vals = np.round(data[valid_mask] / q).astype(np.int32)
    vals = np.clip(vals, -32767, 32767).astype(np.int16)
    out[valid_mask] = vals
    comp = zlib.compress(out.tobytes(order="C"), level=6)
    meta = {"dtype": "i16", "q_m": q, "nodata": int(nodata_i16)}
    return comp, meta

def _pack_grid_f32(data: np.ndarray, valid_mask: np.ndarray) -> Tuple[bytes, Dict[str, Any]]:
    out = data.astype(np.float32, copy=True)
    out[~valid_mask] = np.nan
    comp = zlib.compress(out.tobytes(order="C"), level=6)
    meta = {"dtype": "f32", "nodata": "nan"}
    return comp, meta

def make_bathy_payload_feature_embedded(
    src: rasterio.DatasetReader,
    tb_merc: Tuple[float, float, float, float],
    mode: str,
    grid_size: int,
    quantization_m: float,
    extent: int,
) -> Optional[Dict[str, Any]]:
    mode = (mode or "none").lower()
    if mode == "none":  # explicit no-op
        return None

    out_shape = (grid_size, grid_size)
    window = rasterio.windows.from_bounds(*tb_merc, transform=src.transform)
    data = src.read(1, window=window, out_shape=out_shape, masked=False, resampling=Resampling.bilinear)

    nodata = src.nodata
    if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
        valid = np.isfinite(data)
    else:
        valid = np.isfinite(data) & (data != nodata)

    # MVT requires a geometry; we store payload as one Point feature at tile center
    cx = (tb_merc[0] + tb_merc[2]) * 0.5
    cy = (tb_merc[1] + tb_merc[3]) * 0.5
    px, py = merc_to_tile_coords(cx, cy, tb_merc, extent)

    props: Dict[str, Any] = {
        "mode": mode,
        "grid_w": int(grid_size),
        "grid_h": int(grid_size),
    }

    if mode == "mvt_grid_i16":
        comp, meta = _pack_grid_i16(data, valid, quantization_m)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")
    elif mode == "mvt_grid_f32":
        comp, meta = _pack_grid_f32(data, valid)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")
    elif mode == "mvt_png":
        if not PIL_OK:
            raise RuntimeError("Pillow not installed but bathy mode mvt_png requires it.")
        if valid.any():
            vmin = float(np.nanmin(data[valid]))
            vmax = float(np.nanmax(data[valid]))
        else:
            vmin, vmax = -1.0, 1.0
        if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
            vmin, vmax = -1.0, 1.0
        norm = (data - vmin) / (vmax - vmin)
        norm = np.clip(norm, 0.0, 1.0)
        img = (norm * 255.0).astype(np.uint8)
        img[~valid] = 0
        im = Image.fromarray(img, mode="L")
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=True)
        png_bytes = buf.getvalue()
        comp = zlib.compress(png_bytes, level=6)
        props.update({"dtype": "png8", "vmin": vmin, "vmax": vmax, "nodata_is_zero": True})
        props["png_b64"] = base64.b64encode(comp).decode("ascii")
    else:
        raise ValueError(f"Unknown embedded bathy mode: {mode}")

    return {"geometry": {"type": "Point", "coordinates": [px, py]}, "properties": props}


# ==========================================
# Bathy payload: separate raster MBTiles PNG
# ==========================================

def encode_png8_tile_from_raster(
    data: np.ndarray,
    valid_mask: np.ndarray,
    vmin: float,
    vmax: float,
) -> bytes:
    """Encode a 2D float array into an 8-bit grayscale PNG."""
    if not PIL_OK:
        raise RuntimeError("Pillow not installed but raster MBTiles PNG output requires it.")
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        vmin, vmax = -1.0, 1.0
    norm = (data - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)
    img = (norm * 255.0).astype(np.uint8)
    img[~valid_mask] = 0  # 0 == nodata
    im = Image.fromarray(img, mode="L")
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()

def write_raster_mbtiles_metadata(conn: sqlite3.Connection, bbox_lonlat, zmin, zmax, meta: Dict[str, str]) -> None:
    md = {
        "name": meta.get("name", "Bathymetry raster tiles"),
        "format": "png",
        "minzoom": str(zmin),
        "maxzoom": str(zmax),
        "bounds": ",".join(map(str, bbox_lonlat)),
    }
    md.update(meta)
    mbtiles_put_metadata(conn, md)


# ==============================
# Per-zoom intermediate raster
# ==============================

def build_mosaic_sources(tile_tifs: List[str]) -> List[rasterio.DatasetReader]:
    return [rasterio.open(p) for p in tile_tifs]

def close_mosaic_sources(srcs: List[rasterio.DatasetReader]) -> None:
    for s in srcs:
        try:
            s.close()
        except Exception:
            pass

def ensure_zoom_intermediate(cfg, tile_tifs: List[str], z: int, bbox_merc: Tuple[float, float, float, float]) -> str:
    os.makedirs(cfg.cache_dir, exist_ok=True)
    z_dir = os.path.join(cfg.cache_dir, f"z{z}")
    os.makedirs(z_dir, exist_ok=True)
    out_tif = os.path.join(z_dir, f"dtm_{cfg.dtm_year}_bbox_z{z}_3857.tif")

    if cfg.cache_intermediate and os.path.exists(out_tif):
        return out_tif

    srcs = build_mosaic_sources(tile_tifs)
    try:
        mosaic, mosaic_transform = rio_merge(srcs)
        data = mosaic[0]
        src_crs = srcs[0].crs
        nodata = srcs[0].nodata
        if nodata is None:
            nodata = np.nan
    finally:
        close_mosaic_sources(srcs)

    minx, miny, maxx, maxy = bbox_merc
    res = zoom_target_resolution_m(z, cfg.tile_pixels, tile_size=256)
    width = max(1, int(math.ceil((maxx - minx) / res)))
    height = max(1, int(math.ceil((maxy - miny) / res)))
    dst_transform = from_bounds(minx, miny, maxx, maxy, width=width, height=height)

    dst = np.full((height, width), nodata, dtype=data.dtype)

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


# =========================
# Tile iteration utilities
# =========================

def tile_iter_for_bbox(bbox_lonlat: Tuple[float, float, float, float], z: int) -> Iterable[mercantile.Tile]:
    minlon, minlat, maxlon, maxlat = bbox_lonlat
    return mercantile.tiles(minlon, minlat, maxlon, maxlat, [z])


# =================
# Config dataclass
# =================

@dataclass
class Config:
    # EMODnet DTM discovery
    dtm_year: int
    csw_url: str

    # Paths
    output_mbtiles: str
    cache_dir: str

    # Area / zoom
    bbox_lonlat: Tuple[float, float, float, float]
    zoom_min: int
    zoom_max: int

    # Vector tile generation
    tile_pixels: int
    extent: int
    contour_layer_name: str
    contour_intervals_m: Dict[str, float]

    # Intermediate raster caching
    cache_intermediate: bool

    # Coastline masking
    coast_enabled: bool
    coast_source: str  # osm_land_polygons | emodnet_product
    osm_land_polygons_url: str
    emodnet_vector_url: str
    emodnet_vector_layer: Optional[str]
    coast_max_features: Optional[int]

    # Land mask per-zoom cache
    landmask_cache_enabled: bool
    landmask_all_touched: bool

    # Bathy payload options
    bathy_enabled: bool
    bathy_output: str  # embedded_mvt | raster_mbtiles | none
    bathy_mode: str    # for embedded_mvt: none|mvt_grid_i16|mvt_grid_f32|mvt_png
    bathy_layer_name: str
    bathy_grid_size: int
    bathy_quantization_m: float

    # Raster MBTiles bathy settings
    bathy_raster_mbtiles: str
    bathy_raster_tile_size: int
    bathy_raster_vmin: float
    bathy_raster_vmax: float
    bathy_raster_metadata: Dict[str, str]

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
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    def req(k: str):
        if k not in raw:
            raise ValueError(f"Missing required config key: {k}")
        return raw[k]

    ob = raw.get("obfuscation", {})
    coast = raw.get("coastline_mask", {})
    lmc = raw.get("land_mask_cache", {})
    bathy = raw.get("bathy_data", {})
    raster = bathy.get("raster_mbtiles", {})

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
        osm_land_polygons_url=str(coast.get("osm_land_polygons_url", "https://osmdata.openstreetmap.de/data/land-polygons-split-4326.zip")),
        emodnet_vector_url=str(coast.get("emodnet_vector_url", "")),
        emodnet_vector_layer=coast.get("emodnet_vector_layer", None),
        coast_max_features=(int(coast["max_features"]) if "max_features" in coast and coast["max_features"] is not None else None),

        landmask_cache_enabled=bool(lmc.get("enabled", True)),
        landmask_all_touched=bool(lmc.get("all_touched", False)),

        bathy_enabled=bool(bathy.get("enabled", False)),
        bathy_output=str(bathy.get("output", "none")),
        bathy_mode=str(bathy.get("mode", "none")),
        bathy_layer_name=str(bathy.get("layer_name", "bathy_data")),
        bathy_grid_size=int(bathy.get("grid_size", 32)),
        bathy_quantization_m=float(bathy.get("quantization_m", 0.1)),

        bathy_raster_mbtiles=str(raster.get("output_mbtiles", "")),
        bathy_raster_tile_size=int(raster.get("grid_size", 256)),
        bathy_raster_vmin=float(raster.get("value_scaling", {}).get("vmin", -6000.0)),
        bathy_raster_vmax=float(raster.get("value_scaling", {}).get("vmax", 0.0)),
        bathy_raster_metadata=dict(raster.get("mbtiles_metadata", {})),

        ob_enabled=bool(ob.get("enabled", True)),
        ob_seed=int(ob.get("seed", 12345)),
        ob_jitter_m=dict(ob.get("jitter_m_by_zoom", {})),
        ob_simplify_m=dict(ob.get("simplify_m_by_zoom", {})),
        ob_drop_p=dict(ob.get("drop_probability_by_zoom", {})),
        ob_level_noise_m=dict(ob.get("level_noise_m_by_zoom", {})),

        metadata=dict(raw.get("mbtiles_metadata", {})),
    )


# ==========================
# Main processing pipeline
# ==========================

def run(cfg: Config) -> None:
    rng = random.Random(cfg.ob_seed)

    bbox_merc = lonlat_bbox_to_merc(cfg.bbox_lonlat)
    bbox_geom_merc = box(*bbox_merc)

    # 1) Discover DTM tile URLs
    print("Discovering EMODnet DTM tiles via CSW...")
    tile_urls = csw_search_tiles(cfg.csw_url, cfg.bbox_lonlat, cfg.dtm_year)
    if not tile_urls:
        raise RuntimeError("No EMODnet DTM tile URLs found. Check bbox/year/csw_url.")
    print(f"Found {len(tile_urls)} tile downloads.")

    # 2) Download & extract
    tiles_dir = os.path.join(cfg.cache_dir, "tiles", f"dtm_{cfg.dtm_year}")
    os.makedirs(tiles_dir, exist_ok=True)
    tile_tifs: List[str] = []
    for u in tile_urls:
        print(f"Downloading/extracting: {u}")
        tile_tifs.append(download_and_extract_zip(u, tiles_dir))

    # 3) Load coastline polygons and build spatial index (EPSG:4326)
    land_tree: Optional[STRtree] = None
    land_geoms: List[Polygon | MultiPolygon] = []
    if cfg.coast_enabled:
        if cfg.coast_source == "emodnet_product" and not cfg.emodnet_vector_url:
            raise ValueError("coastline_mask.source=emodnet_product requires coastline_mask.emodnet_vector_url")
        print(f"Loading coastline mask source: {cfg.coast_source}")
        land_tree, land_geoms = load_land_mask_source(cfg, cfg.bbox_lonlat)
        print(f"Loaded land polygons: {len(land_geoms)}")


    # 4) Open vector MBTiles + metadata
    vconn = mbtiles_open(cfg.output_mbtiles)
    vmd = {
        "name": cfg.metadata.get("name", "Bathymetry contours (EMODnet-derived)"),
        "format": "pbf",
        "minzoom": str(cfg.zoom_min),
        "maxzoom": str(cfg.zoom_max),
        "bounds": ",".join(map(str, cfg.bbox_lonlat)),
        "attribution": cfg.metadata.get("attribution", "Derived from EMODnet Bathymetry DTM (CC BY 4.0)."),
    }
    vmd.update(cfg.metadata)
    mbtiles_put_metadata(vconn, vmd)

    # 5) Optional raster MBTiles for bathy field
    rconn: Optional[sqlite3.Connection] = None
    if cfg.bathy_enabled and cfg.bathy_output == "raster_mbtiles":
        if not cfg.bathy_raster_mbtiles:
            raise ValueError("bathy_data.output=raster_mbtiles requires bathy_data.raster_mbtiles.output_mbtiles")
        if not PIL_OK:
            raise RuntimeError("Pillow not installed but raster_mbtiles output requires it.")
        rconn = mbtiles_open(cfg.bathy_raster_mbtiles)
        rmd = {
            "name": cfg.bathy_raster_metadata.get("name", "Bathymetry raster tiles"),
            "description": cfg.bathy_raster_metadata.get("description", "Gridded bathymetry tiles encoded as PNG8"),
            "format": "png",
            "minzoom": str(cfg.zoom_min),
            "maxzoom": str(cfg.zoom_max),
            "bounds": ",".join(map(str, cfg.bbox_lonlat)),
            "bathy_encoding": "png8",
            "bathy_vmin": str(cfg.bathy_raster_vmin),
            "bathy_vmax": str(cfg.bathy_raster_vmax),
            "bathy_tile_size": str(cfg.bathy_raster_tile_size),
        }
        rmd.update(cfg.bathy_raster_metadata)
        mbtiles_put_metadata(rconn, rmd)

    total_tiles = 0
    written_vec = 0
    written_ras = 0
    t0 = time.time()

    for z in range(cfg.zoom_min, cfg.zoom_max + 1):
        interval = _value_for_zoom(cfg.contour_intervals_m, z, None)
        if interval is None or interval <= 0:
            print(f"[z={z}] skip (no contour interval configured)")
            continue

        # Per-zoom intermediate raster
        inter_tif = ensure_zoom_intermediate(cfg, tile_tifs, z, bbox_merc)

        # Per-zoom landmask cache (optional)
        landmask_tif = None
        if cfg.coast_enabled and land_tree is not None:
            landmask_tif = ensure_landmask_zoom(cfg, z, inter_tif, bbox_merc, land_tree, land_geoms)

        # Obfuscation params per zoom
        sigma_jitter = float(_value_for_zoom(cfg.ob_jitter_m, z, 0.0) or 0.0)
        tol_simplify = float(_value_for_zoom(cfg.ob_simplify_m, z, 0.0) or 0.0)
        drop_p = float(_value_for_zoom(cfg.ob_drop_p, z, 0.0) or 0.0)
        level_noise = float(_value_for_zoom(cfg.ob_level_noise_m, z, 0.0) or 0.0)

        print(
            f"[z={z}] interval={interval}m intermediate={os.path.basename(inter_tif)} "
            f"landmask={'ON' if landmask_tif else 'OFF'} "
            f"obfuscation={'ON' if cfg.ob_enabled else 'OFF'} "
            f"bathy_output={cfg.bathy_output}"
        )

        tiles = list(tile_iter_for_bbox(cfg.bbox_lonlat, z))
        total_tiles += len(tiles)

        with rasterio.open(inter_tif) as src, (rasterio.open(landmask_tif) if landmask_tif else nullcontext()) as msrc:  # type: ignore
            nodata = src.nodata

            for t in tiles:
                x, y = t.x, t.y
                tb = tile_bounds_merc(z, x, y)
                tile_geom = box(*tb)
                if not tile_geom.intersects(bbox_geom_merc):
                    continue

                # -------------
                # Contour input
                # -------------
                out_shape = (cfg.tile_pixels, cfg.tile_pixels)
                window = rasterio.windows.from_bounds(*tb, transform=src.transform)
                try:
                    data = src.read(1, window=window, out_shape=out_shape, masked=False, resampling=Resampling.bilinear)
                except Exception:
                    continue

                if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                    nodata_mask = ~np.isfinite(data)
                else:
                    nodata_mask = (data == nodata) | ~np.isfinite(data)

                # Apply landmask if available (window read is cheap)
                if landmask_tif and msrc is not None:
                    mw = rasterio.windows.from_bounds(*tb, transform=msrc.transform)
                    land = msrc.read(1, window=mw, out_shape=out_shape, masked=False, resampling=Resampling.nearest)
                    nodata_mask = nodata_mask | (land.astype(np.uint8) == 1)

                levels = compute_contour_levels(data, float(interval), nodata_mask)
                raw_contours = contours_from_array(data, tb, levels, nodata_mask) if levels else []

                contour_features: List[Dict[str, Any]] = []
                for lev, ls in raw_contours:
                    ls2 = clip_by_rect(ls, tb[0], tb[1], tb[2], tb[3])
                    if ls2.is_empty:
                        continue
                    lines = [ls2] if isinstance(ls2, LineString) else list(ls2.geoms)

                    for line in lines:
                        if line.length <= 0:
                            continue

                        if cfg.ob_enabled:
                            if drop_p > 0 and rng.random() < drop_p:
                                continue
                            line = add_jitter_linestring(line, sigma_jitter, rng)
                            line = simplify_linestring(line, tol_simplify)
                            if line.is_empty or len(line.coords) < 2:
                                continue

                        lev2 = float(lev)
                        if cfg.ob_enabled and level_noise > 0:
                            lev2 += rng.gauss(0.0, level_noise)

                        geom = {
                            "type": "LineString",
                            "coordinates": [merc_to_tile_coords(xx, yy, tb, cfg.extent) for (xx, yy) in line.coords],
                        }
                        contour_features.append({
                            "geometry": geom,
                            "properties": {"depth_m": round(lev2, 3), "interval_m": float(interval), "z": int(z)},
                        })

                # -------------------------
                # Bathy payload (two paths)
                # -------------------------
                bathy_feature = None
                if cfg.bathy_enabled and cfg.bathy_output == "embedded_mvt":
                    try:
                        bathy_feature = make_bathy_payload_feature_embedded(
                            src=src,
                            tb_merc=tb,
                            mode=cfg.bathy_mode,
                            grid_size=cfg.bathy_grid_size,
                            quantization_m=cfg.bathy_quantization_m,
                            extent=cfg.extent,
                        )
                    except Exception:
                        bathy_feature = None

                # Separate raster MBTiles tile
                if rconn is not None:
                    # Read at raster tile size, apply same landmask rule (set nodata pixels to 0)
                    rshape = (cfg.bathy_raster_tile_size, cfg.bathy_raster_tile_size)
                    rw = rasterio.windows.from_bounds(*tb, transform=src.transform)
                    try:
                        rdata = src.read(1, window=rw, out_shape=rshape, masked=False, resampling=Resampling.bilinear)
                    except Exception:
                        rdata = None
                    if rdata is not None:
                        if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                            valid = np.isfinite(rdata)
                        else:
                            valid = np.isfinite(rdata) & (rdata != nodata)

                        # Landmask for raster tile (nearest)
                        if landmask_tif and msrc is not None:
                            rmw = rasterio.windows.from_bounds(*tb, transform=msrc.transform)
                            rland = msrc.read(1, window=rmw, out_shape=rshape, masked=False, resampling=Resampling.nearest)
                            valid = valid & (rland.astype(np.uint8) == 0)

                        png = encode_png8_tile_from_raster(rdata, valid, cfg.bathy_raster_vmin, cfg.bathy_raster_vmax)
                        mbtiles_put_tile(rconn, z, x, tms_y(z, y), png)
                        written_ras += 1

                # If no content, skip writing vector tile
                if not contour_features and bathy_feature is None:
                    continue

                layers: Dict[str, Dict[str, Any]] = {}
                if contour_features:
                    layers[cfg.contour_layer_name] = {"features": contour_features, "extent": cfg.extent, "version": 2}
                if bathy_feature is not None:
                    layers[cfg.bathy_layer_name] = {"features": [bathy_feature], "extent": cfg.extent, "version": 2}

                tile_blob = encode_mvt(layers)
                mbtiles_put_tile(vconn, z, x, tms_y(z, y), tile_blob)
                written_vec += 1

        mbtiles_commit(vconn)
        if rconn is not None:
            mbtiles_commit(rconn)

    dt = time.time() - t0
    print(f"Done. total_tiles={total_tiles}, vector_tiles={written_vec}, raster_tiles={written_ras}, elapsed={dt:.1f}s")

    vconn.commit()
    vconn.close()
    if rconn is not None:
        rconn.commit()
        rconn.close()


# ==========================
# Small contextmanager helper
# ==========================

class nullcontext:
    def __enter__(self):  # noqa
        return None
    def __exit__(self, exc_type, exc, tb):  # noqa
        return False


def main(argv: List[str]) -> int:
    if len(argv) != 2:
        print("Usage: python obfuscate_bathy_tiles.py <config.json>", file=sys.stderr)
        return 2
    cfg = load_config(argv[1])
    run(cfg)
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
