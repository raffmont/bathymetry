#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
obfuscate_bathy_tiles.py

Vector + Raster tile generator for EMODnet bathymetry.

Outputs
-------
1) Vector MBTiles (MVT/PBF): bathymetric contours (+ optional embedded bathy payload)
2) Optional Raster MBTiles (PNG): bathymetry field (png8)
3) Optional Raster MBTiles (PNG): derived rasters (hillshade, slope)
4) TileJSON files for produced MBTiles
5) (Separately) demo HTML files (MapLibre + OpenLayers) can consume TileJSON.

Key features (requested)
------------------------
- CLI overrides (zoom/bbox/disable obfuscation/output paths/etc.)
- Bathy payload either embedded in MVT OR separate raster MBTiles OR none
- Coastline mask source: OSM land polygons OR EMODnet coastline vector product (direct URL)
- Per-zoom landmask GeoTIFF cache aligned to per-zoom intermediate raster
- Optional derived rasters: hillshade and slope (written as raster MBTiles PNG8)
- TileJSON generation for vector + raster MBTiles

Usage
-----
python obfuscate_bathy_tiles.py config.json [--zoom-min 7 --zoom-max 14 --no-obfuscation ...]

Notes
-----
- This script is designed to be "loadable" and reasonably robust, but EMODnet endpoints/products
  can change. If CSW search or a coastline URL changes, adjust the config accordingly.
"""

from __future__ import annotations

import argparse
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
from shapely.geometry import LineString, Polygon, MultiPolygon, box, shape
from shapely.ops import clip_by_rect
from shapely.strtree import STRtree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mapbox_vector_tile

try:
    from PIL import Image
    PIL_OK = True
except Exception:
    PIL_OK = False


# -----------------------------
# Helpers: zoom-range dictionaries
# -----------------------------
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


# -----------------------------
# MBTiles helpers
# -----------------------------
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


# -----------------------------
# WebMercator helpers
# -----------------------------
WEBMERC_MAX = 20037508.342789244

def lonlat_bbox_to_merc(b: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    minlon, minlat, maxlon, maxlat = b
    tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    minx, miny = tf.transform(minlon, minlat)
    maxx, maxy = tf.transform(maxlon, maxlat)
    return (minx, miny, maxx, maxy)

def tile_bounds_merc(z: int, x: int, y: int) -> Tuple[float, float, float, float]:
    b = mercantile.bounds(x, y, z)
    return lonlat_bbox_to_merc((b.west, b.south, b.east, b.north))

def meters_per_pixel_equator(z: int, tile_size: int = 256) -> float:
    world = 2 * WEBMERC_MAX
    return world / (tile_size * 2**z)

def zoom_target_resolution_m(z: int, tile_pixels: int, tile_size: int = 256) -> float:
    return meters_per_pixel_equator(z, tile_size=tile_size) * (tile_size / tile_pixels)


# -----------------------------
# EMODnet CSW discovery
# -----------------------------
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
          <ogc:PropertyIsLike wildCard="*" singleChar="?" escapeChar="\">
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

def csw_search_tiles(csw_url: str, bbox_lonlat: Tuple[float, float, float, float], year: int,
                     timeout_s: int = 60, max_pages: int = 50) -> List[str]:
    minx, miny, maxx, maxy = bbox_lonlat
    seen, out = set(), []
    url_re = re.compile(r"https://downloads\.emodnet-bathymetry\.eu/[^\s\"<]+?_{}\.tif\.zip".format(year))

    start, max_records = 1, 50
    headers = {"Content-Type": "application/xml"}
    for _ in range(max_pages):
        body = CSW_GETRECORDS_TEMPLATE.format(start=start, max_records=max_records, year=year,
                                             minx=minx, miny=miny, maxx=maxx, maxy=maxy)
        r = requests.post(csw_url, data=body.encode("utf-8"), headers=headers, timeout=timeout_s)
        r.raise_for_status()
        xml = r.text
        found = url_re.findall(xml)
        for u in found:
            if u not in seen:
                seen.add(u); out.append(u)
        if 'nextRecord="0"' in xml or 'numberOfRecordsReturned="0"' in xml:
            break
        if not found and start > 1:
            break
        start += max_records
    return out


# -----------------------------
# Download helpers
# -----------------------------
def download_file(url: str, out_path: str, timeout_s: int = 300) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    if os.path.exists(out_path):
        return
    with requests.get(url, stream=True, timeout=timeout_s) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

def download_and_extract_zip(url: str, out_dir: str, timeout_s: int = 300) -> str:
    os.makedirs(out_dir, exist_ok=True)
    local_zip = os.path.join(out_dir, os.path.basename(url.split("?", 1)[0]))
    download_file(url, local_zip, timeout_s=timeout_s)
    with zipfile.ZipFile(local_zip, "r") as z:
        tifs = [m for m in z.namelist() if m.lower().endswith((".tif", ".tiff"))]
        if not tifs:
            raise RuntimeError(f"No GeoTIFF in {local_zip}")
        member = tifs[0]
        out_tif = os.path.join(out_dir, os.path.basename(member))
        if not os.path.exists(out_tif):
            z.extract(member, out_dir)
            extracted = os.path.join(out_dir, member)
            if extracted != out_tif:
                os.makedirs(os.path.dirname(out_tif), exist_ok=True)
                shutil.move(extracted, out_tif)
        return out_tif

def download_and_prepare_vector_dataset(url: str, cache_dir: str) -> str:
    os.makedirs(cache_dir, exist_ok=True)
    fname = os.path.basename(url.split("?", 1)[0])
    local = os.path.join(cache_dir, fname)
    download_file(url, local, timeout_s=600)
    if fname.lower().endswith(".zip"):
        extract_dir = os.path.join(cache_dir, fname[:-4])
        marker = os.path.join(extract_dir, ".extracted")
        if not os.path.exists(marker):
            os.makedirs(extract_dir, exist_ok=True)
            with zipfile.ZipFile(local, "r") as z:
                z.extractall(extract_dir)
            with open(marker, "w", encoding="utf-8") as f:
                f.write("ok\n")
        return extract_dir
    return local

def find_vector_files(root_dir: str) -> List[str]:
    exts = (".gpkg", ".shp", ".geojson", ".json")
    out = []
    for base, _, files in os.walk(root_dir):
        for fn in files:
            if fn.lower().endswith(exts):
                out.append(os.path.join(base, fn))
    out.sort(key=lambda p: (0 if p.lower().endswith(".gpkg") else 1 if p.lower().endswith(".shp") else 2, p))
    return out


# -----------------------------
# Coastline mask loading
# -----------------------------
def bbox_intersects(a: Tuple[float,float,float,float], b: Tuple[float,float,float,float]) -> bool:
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])

def ensure_osm_land_polygons(cache_dir: str, url: str) -> str:
    root = os.path.join(cache_dir, "coastline_mask", "osm_land_polygons")
    os.makedirs(root, exist_ok=True)
    zip_path = os.path.join(root, os.path.basename(url.split("?", 1)[0]))
    download_file(url, zip_path, timeout_s=600)
    marker = os.path.join(root, ".extracted")
    if not os.path.exists(marker):
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(root)
        with open(marker, "w", encoding="utf-8") as f:
            f.write("ok\n")
    return root

def _transform_polygonal(geom, tf: Transformer):
    if geom.is_empty:
        return geom
    if isinstance(geom, Polygon):
        ext = [tf.transform(x, y) for (x, y) in geom.exterior.coords]
        holes = [[tf.transform(x, y) for (x, y) in r.coords] for r in geom.interiors]
        return Polygon(ext, holes)
    if isinstance(geom, MultiPolygon):
        polys = []
        for p in geom.geoms:
            ext = [tf.transform(x, y) for (x, y) in p.exterior.coords]
            holes = [[tf.transform(x, y) for (x, y) in r.coords] for r in p.interiors]
            polys.append(Polygon(ext, holes))
        return MultiPolygon(polys)
    return geom

def load_polygons_epsg4326(path: str, bbox_lonlat: Tuple[float,float,float,float], layer: Optional[str], max_features: Optional[int]) -> List[Polygon|MultiPolygon]:
    out: List[Polygon|MultiPolygon] = []
    with fiona.open(path, layer=layer) as src:
        src_crs = CRS.from_user_input(src.crs_wkt or src.crs or "EPSG:4326")
        tf = Transformer.from_crs(src_crs, "EPSG:4326", always_xy=True)
        cnt = 0
        for feat in src:
            g = feat.get("geometry")
            if not g:
                continue
            s = shape(g)
            if not isinstance(s, (Polygon, MultiPolygon)) or s.is_empty:
                continue
            if src_crs.to_string() != "EPSG:4326":
                s = _transform_polygonal(s, tf)
            if not bbox_intersects(s.bounds, bbox_lonlat):
                continue
            out.append(s)
            cnt += 1
            if max_features and cnt >= max_features:
                break
    return out

def load_land_polygons(cfg, bbox_lonlat: Tuple[float,float,float,float]) -> Tuple[STRtree, List[Polygon|MultiPolygon]]:
    if not cfg.coast_enabled:
        return STRtree([]), []
    if cfg.coast_source == "osm_land_polygons":
        root = ensure_osm_land_polygons(cfg.cache_dir, cfg.osm_land_polygons_url)
        shps = []
        for base, _, files in os.walk(root):
            for fn in files:
                if fn.lower().endswith(".shp"):
                    shps.append(os.path.join(base, fn))
        geoms: List[Polygon|MultiPolygon] = []
        for shp in shps:
            try:
                with fiona.open(shp) as s:
                    if not bbox_intersects(s.bounds, bbox_lonlat):
                        continue
            except Exception:
                continue
            geoms.extend(load_polygons_epsg4326(shp, bbox_lonlat, None, cfg.coast_max_features))
            if cfg.coast_max_features and len(geoms) >= cfg.coast_max_features:
                geoms = geoms[:cfg.coast_max_features]
                break
        return (STRtree(geoms) if geoms else STRtree([]), geoms)

    if cfg.coast_source == "emodnet_product":
        if not cfg.emodnet_vector_url:
            raise ValueError("coastline_mask.source=emodnet_product requires coastline_mask.emodnet_vector_url")
        vcache = os.path.join(cfg.cache_dir, "coastline_mask", "emodnet_product")
        ds = download_and_prepare_vector_dataset(cfg.emodnet_vector_url, vcache)
        if os.path.isdir(ds):
            cands = find_vector_files(ds)
            if not cands:
                raise RuntimeError(f"No vector files in extracted dataset: {ds}")
            ds = cands[0]
        geoms = load_polygons_epsg4326(ds, bbox_lonlat, cfg.emodnet_vector_layer, cfg.coast_max_features)
        return (STRtree(geoms) if geoms else STRtree([]), geoms)

    raise ValueError(f"Unsupported coastline_mask.source: {cfg.coast_source}")


# -----------------------------
# Per-zoom intermediate raster + landmask cache
# -----------------------------
def build_mosaic_sources(tile_tifs: List[str]) -> List[rasterio.DatasetReader]:
    return [rasterio.open(p) for p in tile_tifs]

def close_mosaic_sources(srcs: List[rasterio.DatasetReader]) -> None:
    for s in srcs:
        try:
            s.close()
        except Exception:
            pass

def ensure_zoom_intermediate(cfg, tile_tifs: List[str], z: int, bbox_merc: Tuple[float,float,float,float]) -> str:
    zdir = os.path.join(cfg.cache_dir, f"z{z}")
    os.makedirs(zdir, exist_ok=True)
    out_tif = os.path.join(zdir, f"dtm_{cfg.dtm_year}_bbox_z{z}_3857.tif")
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
        source=data, destination=dst,
        src_transform=mosaic_transform, src_crs=src_crs,
        dst_transform=dst_transform, dst_crs="EPSG:3857",
        resampling=Resampling.bilinear, src_nodata=nodata, dst_nodata=nodata
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

def ensure_landmask_zoom(cfg, z: int, inter_tif: str, bbox_merc: Tuple[float,float,float,float],
                        land_tree: STRtree, land_geoms: List[Polygon|MultiPolygon]) -> Optional[str]:
    if not (cfg.landmask_cache_enabled and cfg.coast_enabled):
        return None
    zdir = os.path.join(cfg.cache_dir, f"z{z}")
    os.makedirs(zdir, exist_ok=True)
    out_mask = os.path.join(zdir, f"landmask_z{z}_3857.tif")
    if os.path.exists(out_mask):
        return out_mask

    with rasterio.open(inter_tif) as src:
        profile = src.profile.copy()
        profile.update({"count": 1, "dtype": "uint8", "nodata": 0, "compress": "DEFLATE", "predictor": 2})
        H, W = src.height, src.width
        transform = src.transform

        if not land_geoms:
            mask = np.zeros((H, W), dtype=np.uint8)
        else:
            bbox_lonlat = cfg.bbox_lonlat
            candidates = land_tree.query(box(*bbox_lonlat))
            tf = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            minx, miny, maxx, maxy = bbox_merc
            shapes_3857 = []
            for g in candidates:
                if not bbox_intersects(g.bounds, bbox_lonlat):
                    continue
                gg = _transform_polygonal(g, tf)
                gg = clip_by_rect(gg, minx, miny, maxx, maxy)
                if gg.is_empty:
                    continue
                shapes_3857.append(gg)

            mask = rasterize(
                [(g, 1) for g in shapes_3857],
                out_shape=(H, W),
                transform=transform,
                fill=0,
                dtype=np.uint8,
                all_touched=bool(cfg.landmask_all_touched),
            )

        with rasterio.open(out_mask, "w", **profile) as dst:
            dst.write(mask, 1)

    return out_mask


# -----------------------------
# Contours + MVT encoding
# -----------------------------
def compute_contour_levels(data: np.ndarray, interval_m: float, nodata_mask: np.ndarray) -> List[float]:
    valid = data[~nodata_mask]
    if valid.size == 0:
        return []
    vmin, vmax = float(np.nanmin(valid)), float(np.nanmax(valid))
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        return []
    lo = math.floor(vmin / interval_m) * interval_m
    hi = math.ceil(vmax / interval_m) * interval_m
    n = int(round((hi - lo) / interval_m)) + 1
    if n <= 1 or n > 20000:
        return []
    return [lo + i * interval_m for i in range(n)]

def contours_from_array(grid: np.ndarray, extent_merc: Tuple[float,float,float,float], levels: List[float],
                        nodata_mask: np.ndarray) -> List[Tuple[float, LineString]]:
    if not levels:
        return []
    g = grid.copy()
    g[nodata_mask] = np.nan
    minx, miny, maxx, maxy = extent_merc
    h, w = g.shape
    xs = np.linspace(minx, maxx, w)
    ys = np.linspace(miny, maxy, h)
    fig = plt.figure(figsize=(1, 1), dpi=72)
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    try:
        cs = ax.contour(xs, ys, g, levels=levels, linewidths=1)
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

def merc_to_tile_coords(x: float, y: float, bounds: Tuple[float,float,float,float], extent: int) -> Tuple[int, int]:
    minx, miny, maxx, maxy = bounds
    tx = (x - minx) / (maxx - minx)
    ty = (maxy - y) / (maxy - miny)
    px = int(round(tx * extent))
    py = int(round(ty * extent))
    return (max(0, min(extent, px)), max(0, min(extent, py)))

def encode_mvt(layers: Dict[str, Dict[str, Any]]) -> bytes:
    payload = {"layers": []}
    for name, layer in layers.items():
        payload["layers"].append({
            "name": name,
            "features": layer["features"],
            "extent": layer["extent"],
            "version": layer.get("version", 2),
        })
    return mapbox_vector_tile.encode(payload)


# -----------------------------
# Bathy payload: embedded MVT
# -----------------------------
def _pack_grid_i16(data: np.ndarray, valid_mask: np.ndarray, q_m: float) -> Tuple[bytes, Dict[str, Any]]:
    q = float(q_m) if q_m and q_m > 0 else 1.0
    nod = np.int16(-32768)
    out = np.full(data.shape, nod, dtype=np.int16)
    vals = np.round(data[valid_mask] / q).astype(np.int32)
    vals = np.clip(vals, -32767, 32767).astype(np.int16)
    out[valid_mask] = vals
    comp = zlib.compress(out.tobytes(order="C"), level=6)
    return comp, {"dtype": "i16", "q_m": q, "nodata": int(nod)}

def _pack_grid_f32(data: np.ndarray, valid_mask: np.ndarray) -> Tuple[bytes, Dict[str, Any]]:
    out = data.astype(np.float32, copy=True)
    out[~valid_mask] = np.nan
    comp = zlib.compress(out.tobytes(order="C"), level=6)
    return comp, {"dtype": "f32", "nodata": "nan"}

def make_bathy_payload_feature_embedded(src: rasterio.DatasetReader, tb_merc: Tuple[float,float,float,float],
                                       mode: str, grid_size: int, quantization_m: float, extent: int) -> Optional[Dict[str, Any]]:
    mode = (mode or "none").lower()
    if mode == "none":
        return None
    window = rasterio.windows.from_bounds(*tb_merc, transform=src.transform)
    data = src.read(1, window=window, out_shape=(grid_size, grid_size), masked=False, resampling=Resampling.bilinear)
    nodata = src.nodata
    valid = np.isfinite(data) if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)) else (np.isfinite(data) & (data != nodata))

    cx = (tb_merc[0] + tb_merc[2]) * 0.5
    cy = (tb_merc[1] + tb_merc[3]) * 0.5
    px, py = merc_to_tile_coords(cx, cy, tb_merc, extent)

    props: Dict[str, Any] = {"mode": mode, "grid_w": int(grid_size), "grid_h": int(grid_size)}
    if mode == "mvt_grid_i16":
        comp, meta = _pack_grid_i16(data, valid, quantization_m)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")
    elif mode == "mvt_grid_f32":
        comp, meta = _pack_grid_f32(data, valid)
        props.update(meta)
        props["grid_b64"] = base64.b64encode(comp).decode("ascii")
    else:
        raise ValueError(f"Unknown embedded bathy mode: {mode}")

    return {"geometry": {"type": "Point", "coordinates": [px, py]}, "properties": props}


# -----------------------------
# Raster tiles: bathy + derived rasters (PNG8)
# -----------------------------
def encode_png8(data: np.ndarray, valid: np.ndarray, vmin: float, vmax: float) -> bytes:
    if not PIL_OK:
        raise RuntimeError("Pillow required for raster MBTiles outputs")
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        vmin, vmax = -1.0, 1.0
    norm = np.clip((data - vmin) / (vmax - vmin), 0.0, 1.0)
    img = (norm * 255.0).astype(np.uint8)
    img[~valid] = 0
    im = Image.fromarray(img, mode="L")
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()

def hillshade_and_slope(z_arr: np.ndarray, valid: np.ndarray, pixel_size_m: float,
                        azimuth_deg: float, altitude_deg: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    zz = z_arr.astype(np.float32, copy=True)
    zz[~valid] = np.nan
    gy, gx = np.gradient(zz, pixel_size_m, pixel_size_m)

    slope = np.arctan(np.sqrt(gx * gx + gy * gy))
    slope_deg = np.degrees(slope)

    aspect = np.arctan2(-gy, gx)
    az = np.radians(azimuth_deg)
    alt = np.radians(altitude_deg)

    hs = (np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect))
    hs = np.clip(hs, 0.0, 1.0)
    hs255 = (hs * 255.0).astype(np.float32)

    out_valid = np.isfinite(slope_deg) & valid
    slope_deg = np.clip(slope_deg, 0.0, 90.0).astype(np.float32)
    return hs255, slope_deg, out_valid


# -----------------------------
# TileJSON generation
# -----------------------------
def write_tilejson(path: str, tilejson: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(tilejson, f, indent=2)

def make_tilejson(name: str, fmt: str, bbox_lonlat, zmin, zmax, tiles_url_template: str, extra: Dict[str, Any]) -> Dict[str, Any]:
    tj = {
        "tilejson": "3.0.0",
        "name": name,
        "format": fmt,
        "scheme": "xyz",
        "minzoom": int(zmin),
        "maxzoom": int(zmax),
        "bounds": list(map(float, bbox_lonlat)),
        "tiles": [tiles_url_template],
    }
    tj.update(extra or {})
    return tj


# -----------------------------
# Config model
# -----------------------------
@dataclass
class Config:
    dtm_year: int
    csw_url: str

    output_mbtiles: str
    cache_dir: str

    bbox_lonlat: Tuple[float, float, float, float]
    zoom_min: int
    zoom_max: int

    tile_pixels: int
    extent: int
    contour_layer_name: str
    contour_intervals_m: Dict[str, float]

    cache_intermediate: bool

    coast_enabled: bool
    coast_source: str
    osm_land_polygons_url: str
    emodnet_vector_url: str
    emodnet_vector_layer: Optional[str]
    coast_max_features: Optional[int]

    landmask_cache_enabled: bool
    landmask_all_touched: bool

    bathy_enabled: bool
    bathy_output: str
    bathy_mode: str
    bathy_layer_name: str
    bathy_grid_size: int
    bathy_quantization_m: float

    bathy_raster_mbtiles: str
    bathy_raster_tile_size: int
    bathy_raster_vmin: float
    bathy_raster_vmax: float
    bathy_raster_metadata: Dict[str, str]

    derived_enabled: bool
    derived_hillshade_mbtiles: str
    derived_slope_mbtiles: str
    derived_tile_size: int
    derived_hs_azimuth: float
    derived_hs_altitude: float
    derived_metadata: Dict[str, str]

    tilejson_enabled: bool
    tilejson_base_url: str
    tilejson_out_dir: str

    ob_enabled: bool
    ob_seed: int
    ob_jitter_m: Dict[str, float]
    ob_simplify_m: Dict[str, float]
    ob_drop_p: Dict[str, float]
    ob_level_noise_m: Dict[str, float]

    metadata: Dict[str, str]


def load_config(path: str) -> Config:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    def req(k: str):
        if k not in raw:
            raise ValueError(f"Missing required config key: {k}")
        return raw[k]

    coast = raw.get("coastline_mask", {})
    lmc = raw.get("land_mask_cache", {})
    bathy = raw.get("bathy_data", {})
    raster = bathy.get("raster_mbtiles", {})
    derived = raw.get("derived_rasters", {})
    tilejson = raw.get("tilejson", {})
    ob = raw.get("obfuscation", {})

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

        derived_enabled=bool(derived.get("enabled", False)),
        derived_hillshade_mbtiles=str(derived.get("hillshade_mbtiles", "")),
        derived_slope_mbtiles=str(derived.get("slope_mbtiles", "")),
        derived_tile_size=int(derived.get("tile_size", 256)),
        derived_hs_azimuth=float(derived.get("hillshade", {}).get("azimuth_deg", 315.0)),
        derived_hs_altitude=float(derived.get("hillshade", {}).get("altitude_deg", 45.0)),
        derived_metadata=dict(derived.get("mbtiles_metadata", {})),

        tilejson_enabled=bool(tilejson.get("enabled", True)),
        tilejson_base_url=str(tilejson.get("base_url", "http://localhost:8080/tiles")),
        tilejson_out_dir=str(tilejson.get("out_dir", "./data/output/tilejson")),

        ob_enabled=bool(ob.get("enabled", True)),
        ob_seed=int(ob.get("seed", 12345)),
        ob_jitter_m=dict(ob.get("jitter_m_by_zoom", {})),
        ob_simplify_m=dict(ob.get("simplify_m_by_zoom", {})),
        ob_drop_p=dict(ob.get("drop_probability_by_zoom", {})),
        ob_level_noise_m=dict(ob.get("level_noise_m_by_zoom", {})),

        metadata=dict(raw.get("mbtiles_metadata", {})),
    )


# -----------------------------
# CLI overrides
# -----------------------------
def apply_cli_overrides(cfg: Config, args: argparse.Namespace) -> Config:
    d = cfg.__dict__.copy()
    if args.zoom_min is not None:
        d["zoom_min"] = args.zoom_min
    if args.zoom_max is not None:
        d["zoom_max"] = args.zoom_max
    if args.bbox is not None:
        d["bbox_lonlat"] = tuple(map(float, args.bbox))
    if args.no_obfuscation:
        d["ob_enabled"] = False
    if args.vector_mbtiles is not None:
        d["output_mbtiles"] = args.vector_mbtiles
    if args.bathy_output is not None:
        d["bathy_output"] = args.bathy_output
    if args.raster_mbtiles is not None:
        d["bathy_raster_mbtiles"] = args.raster_mbtiles
    if args.coast_source is not None:
        d["coast_source"] = args.coast_source
    if args.emodnet_coast_url is not None:
        d["emodnet_vector_url"] = args.emodnet_coast_url
    return Config(**d)


# -----------------------------
# Tile iteration
# -----------------------------
def tile_iter_for_bbox(bbox_lonlat: Tuple[float,float,float,float], z: int) -> Iterable[mercantile.Tile]:
    minlon, minlat, maxlon, maxlat = bbox_lonlat
    return mercantile.tiles(minlon, minlat, maxlon, maxlat, [z])


# -----------------------------
# Main pipeline
# -----------------------------
def run(cfg: Config) -> None:
    rng = random.Random(cfg.ob_seed)
    bbox_merc = lonlat_bbox_to_merc(cfg.bbox_lonlat)
    bbox_geom_merc = box(*bbox_merc)

    print("Discovering EMODnet DTM tiles via CSW...")
    tile_urls = csw_search_tiles(cfg.csw_url, cfg.bbox_lonlat, cfg.dtm_year)
    if not tile_urls:
        raise RuntimeError("No EMODnet DTM tile URLs found. Check bbox/year/csw_url.")
    print(f"Found {len(tile_urls)} tile downloads.")

    tiles_dir = os.path.join(cfg.cache_dir, "tiles", f"dtm_{cfg.dtm_year}")
    os.makedirs(tiles_dir, exist_ok=True)
    tile_tifs: List[str] = []
    for u in tile_urls:
        print(f"Downloading/extracting: {u}")
        tile_tifs.append(download_and_extract_zip(u, tiles_dir))

    land_tree, land_geoms = load_land_polygons(cfg, cfg.bbox_lonlat) if cfg.coast_enabled else (STRtree([]), [])
    if cfg.coast_enabled:
        print(f"Loaded land polygons: {len(land_geoms)} (source={cfg.coast_source})")

    vconn = mbtiles_open(cfg.output_mbtiles)
    vmd = {
        "name": cfg.metadata.get("name", "Bathymetry contours (EMODnet-derived)"),
        "format": "pbf",
        "minzoom": str(cfg.zoom_min),
        "maxzoom": str(cfg.zoom_max),
        "bounds": ",".join(map(str, cfg.bbox_lonlat)),
        "attribution": cfg.metadata.get("attribution", "Derived from EMODnet Bathymetry DTM (CC BY 4.0)."),
        "description": cfg.metadata.get("description", ""),
    }
    vmd.update(cfg.metadata)
    mbtiles_put_metadata(vconn, vmd)

    rconn = None
    if cfg.bathy_enabled and cfg.bathy_output == "raster_mbtiles":
        if not PIL_OK:
            raise RuntimeError("Pillow not installed but raster_mbtiles output requires it.")
        if not cfg.bathy_raster_mbtiles:
            raise ValueError("bathy_data.output=raster_mbtiles requires bathy_data.raster_mbtiles.output_mbtiles")
        rconn = mbtiles_open(cfg.bathy_raster_mbtiles)
        rmd = {
            "name": cfg.bathy_raster_metadata.get("name", "Bathymetry raster tiles (png8)"),
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

    hsconn = slconn = None
    if cfg.derived_enabled:
        if not PIL_OK:
            raise RuntimeError("Pillow not installed but derived rasters require it.")
        if not cfg.derived_hillshade_mbtiles or not cfg.derived_slope_mbtiles:
            raise ValueError("derived_rasters.enabled=true requires hillshade_mbtiles and slope_mbtiles.")
        hsconn = mbtiles_open(cfg.derived_hillshade_mbtiles)
        slconn = mbtiles_open(cfg.derived_slope_mbtiles)
        dmd = {
            "format": "png",
            "minzoom": str(cfg.zoom_min),
            "maxzoom": str(cfg.zoom_max),
            "bounds": ",".join(map(str, cfg.bbox_lonlat)),
            "tile_size": str(cfg.derived_tile_size),
        }
        md_hs = {"name": "Hillshade (png8)", "description": "Hillshade derived from bathymetry",
                 "hs_azimuth_deg": str(cfg.derived_hs_azimuth), "hs_altitude_deg": str(cfg.derived_hs_altitude), **dmd}
        md_sl = {"name": "Slope (png8)", "description": "Slope degrees derived from bathymetry", "slope_units": "degrees", **dmd}
        md_hs.update(cfg.derived_metadata)
        md_sl.update(cfg.derived_metadata)
        mbtiles_put_metadata(hsconn, md_hs)
        mbtiles_put_metadata(slconn, md_sl)

    total_tiles = written_vec = written_ras = written_hs = written_sl = 0
    t0 = time.time()

    for z in range(cfg.zoom_min, cfg.zoom_max + 1):
        interval = _value_for_zoom(cfg.contour_intervals_m, z, None)
        if interval is None or interval <= 0:
            print(f"[z={z}] skip (no contour interval configured)")
            continue

        inter_tif = ensure_zoom_intermediate(cfg, tile_tifs, z, bbox_merc)
        landmask_tif = ensure_landmask_zoom(cfg, z, inter_tif, bbox_merc, land_tree, land_geoms)

        sigma_jitter = float(_value_for_zoom(cfg.ob_jitter_m, z, 0.0) or 0.0)
        tol_simplify = float(_value_for_zoom(cfg.ob_simplify_m, z, 0.0) or 0.0)
        drop_p = float(_value_for_zoom(cfg.ob_drop_p, z, 0.0) or 0.0)
        level_noise = float(_value_for_zoom(cfg.ob_level_noise_m, z, 0.0) or 0.0)

        print(f"[z={z}] interval={interval}m inter={os.path.basename(inter_tif)} landmask={'ON' if landmask_tif else 'OFF'}")

        tiles = list(tile_iter_for_bbox(cfg.bbox_lonlat, z))
        total_tiles += len(tiles)

        with rasterio.open(inter_tif) as src:
            msrc = rasterio.open(landmask_tif) if landmask_tif else None
            nodata = src.nodata

            for t in tiles:
                x, y = t.x, t.y
                tb = tile_bounds_merc(z, x, y)
                if not box(*tb).intersects(bbox_geom_merc):
                    continue

                out_shape = (cfg.tile_pixels, cfg.tile_pixels)
                w = rasterio.windows.from_bounds(*tb, transform=src.transform)
                try:
                    data = src.read(1, window=w, out_shape=out_shape, masked=False, resampling=Resampling.bilinear)
                except Exception:
                    continue

                if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                    valid = np.isfinite(data)
                else:
                    valid = np.isfinite(data) & (data != nodata)

                if msrc is not None:
                    mw = rasterio.windows.from_bounds(*tb, transform=msrc.transform)
                    land = msrc.read(1, window=mw, out_shape=out_shape, masked=False, resampling=Resampling.nearest)
                    valid = valid & (land.astype(np.uint8) == 0)

                nodata_mask = ~valid

                levels = compute_contour_levels(data, float(interval), nodata_mask)
                raw = contours_from_array(data, tb, levels, nodata_mask) if levels else []

                contour_features: List[Dict[str, Any]] = []
                for lev, ls in raw:
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

                        geom = {"type": "LineString",
                                "coordinates": [merc_to_tile_coords(xx, yy, tb, cfg.extent) for (xx, yy) in line.coords]}
                        contour_features.append({"geometry": geom, "properties": {"depth_m": round(lev2, 3), "interval_m": float(interval), "z": int(z)}})

                bathy_feature = None
                if cfg.bathy_enabled and cfg.bathy_output == "embedded_mvt":
                    try:
                        bathy_feature = make_bathy_payload_feature_embedded(src, tb, cfg.bathy_mode, cfg.bathy_grid_size, cfg.bathy_quantization_m, cfg.extent)
                    except Exception:
                        bathy_feature = None

                if rconn is not None:
                    rshape = (cfg.bathy_raster_tile_size, cfg.bathy_raster_tile_size)
                    rw = rasterio.windows.from_bounds(*tb, transform=src.transform)
                    rdata = src.read(1, window=rw, out_shape=rshape, masked=False, resampling=Resampling.bilinear)
                    if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                        rvalid = np.isfinite(rdata)
                    else:
                        rvalid = np.isfinite(rdata) & (rdata != nodata)
                    if msrc is not None:
                        rmw = rasterio.windows.from_bounds(*tb, transform=msrc.transform)
                        rland = msrc.read(1, window=rmw, out_shape=rshape, masked=False, resampling=Resampling.nearest)
                        rvalid = rvalid & (rland.astype(np.uint8) == 0)
                    png = encode_png8(rdata, rvalid, cfg.bathy_raster_vmin, cfg.bathy_raster_vmax)
                    mbtiles_put_tile(rconn, z, x, tms_y(z, y), png)
                    written_ras += 1

                if hsconn is not None and slconn is not None:
                    dshape = (cfg.derived_tile_size, cfg.derived_tile_size)
                    dw = rasterio.windows.from_bounds(*tb, transform=src.transform)
                    ddata = src.read(1, window=dw, out_shape=dshape, masked=False, resampling=Resampling.bilinear)
                    if nodata is None or (isinstance(nodata, float) and math.isnan(nodata)):
                        dvalid = np.isfinite(ddata)
                    else:
                        dvalid = np.isfinite(ddata) & (ddata != nodata)
                    if msrc is not None:
                        dmw = rasterio.windows.from_bounds(*tb, transform=msrc.transform)
                        dland = msrc.read(1, window=dmw, out_shape=dshape, masked=False, resampling=Resampling.nearest)
                        dvalid = dvalid & (dland.astype(np.uint8) == 0)

                    pixel_size = zoom_target_resolution_m(z, cfg.derived_tile_size, tile_size=256)
                    hs255, slope_deg, dvalid2 = hillshade_and_slope(ddata, dvalid, pixel_size, cfg.derived_hs_azimuth, cfg.derived_hs_altitude)

                    png_hs = encode_png8(hs255, dvalid2, 0.0, 255.0)
                    mbtiles_put_tile(hsconn, z, x, tms_y(z, y), png_hs)
                    written_hs += 1

                    png_sl = encode_png8(slope_deg, dvalid2, 0.0, 90.0)
                    mbtiles_put_tile(slconn, z, x, tms_y(z, y), png_sl)
                    written_sl += 1

                if not contour_features and bathy_feature is None:
                    continue

                layers: Dict[str, Dict[str, Any]] = {}
                if contour_features:
                    layers[cfg.contour_layer_name] = {"features": contour_features, "extent": cfg.extent, "version": 2}
                if bathy_feature is not None:
                    layers[cfg.bathy_layer_name] = {"features": [bathy_feature], "extent": cfg.extent, "version": 2}

                blob = encode_mvt(layers)
                mbtiles_put_tile(vconn, z, x, tms_y(z, y), blob)
                written_vec += 1

            if msrc is not None:
                msrc.close()

        mbtiles_commit(vconn)
        if rconn is not None:
            mbtiles_commit(rconn)
        if hsconn is not None:
            mbtiles_commit(hsconn)
        if slconn is not None:
            mbtiles_commit(slconn)

    dt = time.time() - t0
    print(f"Done. total_tiles={total_tiles}, vector={written_vec}, bathy_png={written_ras}, hillshade={written_hs}, slope={written_sl}, elapsed={dt:.1f}s")

    vconn.close()
    if rconn is not None:
        rconn.close()
    if hsconn is not None:
        hsconn.close()
    if slconn is not None:
        slconn.close()

    if cfg.tilejson_enabled:
        out_dir = cfg.tilejson_out_dir
        os.makedirs(out_dir, exist_ok=True)

        vec_url = f"{cfg.tilejson_base_url}/vector/{{z}}/{{x}}/{{y}}.pbf"
        tj_vec = make_tilejson("bathy_contours", "pbf", cfg.bbox_lonlat, cfg.zoom_min, cfg.zoom_max, vec_url, {
            "vector_layers": [{"id": cfg.contour_layer_name, "description": "Bathymetric contours"}],
        })
        write_tilejson(os.path.join(out_dir, "bathy_contours.tilejson"), tj_vec)

        if cfg.bathy_enabled and cfg.bathy_output == "raster_mbtiles":
            ras_url = f"{cfg.tilejson_base_url}/bathy/{{z}}/{{x}}/{{y}}.png"
            tj_ras = make_tilejson("bathy_raster", "png", cfg.bbox_lonlat, cfg.zoom_min, cfg.zoom_max, ras_url, {
                "bathy_encoding": "png8",
                "bathy_vmin": cfg.bathy_raster_vmin,
                "bathy_vmax": cfg.bathy_raster_vmax,
                "tile_size": cfg.bathy_raster_tile_size,
            })
            write_tilejson(os.path.join(out_dir, "bathy_raster.tilejson"), tj_ras)

        if cfg.derived_enabled:
            hs_url = f"{cfg.tilejson_base_url}/hillshade/{{z}}/{{x}}/{{y}}.png"
            sl_url = f"{cfg.tilejson_base_url}/slope/{{z}}/{{x}}/{{y}}.png"
            tj_hs = make_tilejson("hillshade", "png", cfg.bbox_lonlat, cfg.zoom_min, cfg.zoom_max, hs_url, {"tile_size": cfg.derived_tile_size})
            tj_sl = make_tilejson("slope", "png", cfg.bbox_lonlat, cfg.zoom_min, cfg.zoom_max, sl_url, {"tile_size": cfg.derived_tile_size})
            write_tilejson(os.path.join(out_dir, "hillshade.tilejson"), tj_hs)
            write_tilejson(os.path.join(out_dir, "slope.tilejson"), tj_sl)


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate obfuscated bathy contour vector tiles + optional raster tiles.")
    p.add_argument("config", help="Path to JSON config.")
    p.add_argument("--zoom-min", type=int, default=None)
    p.add_argument("--zoom-max", type=int, default=None)
    p.add_argument("--bbox", type=float, nargs=4, default=None, metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"))
    p.add_argument("--no-obfuscation", action="store_true")
    p.add_argument("--vector-mbtiles", default=None, help="Override output vector MBTiles path")
    p.add_argument("--bathy-output", default=None, choices=["none", "embedded_mvt", "raster_mbtiles"], help="Override bathy output strategy")
    p.add_argument("--raster-mbtiles", default=None, help="Override bathy raster MBTiles path")
    p.add_argument("--coast-source", default=None, choices=["osm_land_polygons", "emodnet_product"], help="Override coastline mask source")
    p.add_argument("--emodnet-coast-url", default=None, help="Override EMODnet coastline product URL")
    return p.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv[1:])
    cfg = load_config(args.config)
    cfg = apply_cli_overrides(cfg, args)
    run(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
