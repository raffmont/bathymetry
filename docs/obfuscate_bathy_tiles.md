# obfuscate_bathy_tiles

Generate (optionally obfuscated) bathymetric contour vector tiles (MVT) in an MBTiles container, starting from EMODnet Bathymetry DTM tiles.

This version supports:
- Auto-discovery + download of EMODnet DTM tiles for a bbox (via EMODnet GeoNetwork CSW)
- Mosaic on the fly + per-zoom intermediate raster cache (EPSG:3857)
- Contour intervals configurable per zoom or zoom-range
- Optional obfuscation (jitter/simplify/drop + depth noise)
- Land/coastline masking (OSM land polygons **or** EMODnet coastline products)
- **Per-zoom land mask raster cache** (so masking is fast at tile time)
- Optional per-tile bathymetry payload in multiple modes:
  - Embedded inside MVT
  - **OR** written as a **separate raster MBTiles** (PNG tiles)

> **Licensing / attribution**  
> EMODnet data products are generally CC BY 4.0; you are responsible for respecting EMODnet terms and proper attribution.

---

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Run

```bash
python obfuscate_bathy_tiles.py config.json
```

Outputs:
- `output_mbtiles` (vector MBTiles with MVT tiles)
- Optionally: `bathy_data.raster_mbtiles.output_mbtiles` (raster MBTiles with PNG tiles)

---

## Pipeline overview

### 1) Discover DTM tiles intersecting bbox
Uses EMODnet GeoNetwork CSW search for records matching *EMODnet Digital Bathymetry (DTM YEAR) - Tile* and a bbox filter.

Config keys:
- `dtm_year`
- `csw_url`
- `bbox_lonlat`

### 2) Download + extract tiles (cached)
Tiles are cached under:
- `{cache_dir}/tiles/dtm_{dtm_year}/`

### 3) Mosaic + per-zoom intermediate raster cache
For each zoom `z`, we create or reuse an intermediate raster:
- mosaic all DTM tiles
- warp+clip to EPSG:3857 bbox
- resample to zoom-appropriate resolution
- write GeoTIFF cache:
  - `{cache_dir}/z{z}/dtm_{year}_bbox_z{z}_3857.tif`

Config keys:
- `cache_intermediate` (true/false)
- `tile_pixels` (affects chosen warp resolution)

### 4) Coastline / land masking (two sources)

#### A) OSM land polygons (recommended default)
Downloads and caches the OSM “land polygons split 4326” dataset, loads only polygons intersecting bbox, and uses them to build a per-zoom land mask raster.

```json
"coastline_mask": {
  "enabled": true,
  "source": "osm_land_polygons",
  "osm_land_polygons_url": "https://osmdata.openstreetmap.de/data/land-polygons-split-4326.zip"
}
```

#### B) EMODnet coastline products
Sometimes you may prefer a coastline dataset distributed by EMODnet (e.g., generalized coastlines, regional products, etc.).

Because EMODnet product packaging and endpoints can vary, the tool supports a **direct download URL** you provide to one of:
- a zipped shapefile (ZIP containing .shp/.dbf/.shx/.prj)
- a GeoPackage (`.gpkg`)
- a GeoJSON (`.geojson` / `.json`)
- any Fiona-readable vector dataset

```json
"coastline_mask": {
  "enabled": true,
  "source": "emodnet_product",
  "emodnet_vector_url": "https://.../coastline.zip",
  "emodnet_vector_layer": null
}
```

Notes:
- `emodnet_vector_layer` is optional; only needed for multi-layer containers (e.g., GeoPackage).

### 5) Per-zoom land mask raster cache (fast masking)
Instead of rasterizing polygons for every single tile, the tool builds **one land mask raster per zoom**, aligned to the intermediate raster grid.

Cache file:
- `{cache_dir}/z{z}/landmask_z{z}_3857.tif`

At tile time:
- read the landmask window for the tile
- treat land pixels as `nodata` (so contours aren’t produced on land)

Config keys:
```json
"land_mask_cache": {
  "enabled": true,
  "all_touched": false
}
```

- `all_touched=false` (conservative): a pixel is land only if the polygon covers the pixel center (default rasterize behavior)
- `all_touched=true`: more aggressive masking along coasts

### 6) Contour intervals per zoom or zoom-range
```json
"contour_intervals_m": {
  "7-9": 50,
  "10-11": 20,
  "12-14": 10
}
```

### 7) Optional obfuscation
```json
"obfuscation": {
  "enabled": true,
  "seed": 20260201,
  "jitter_m_by_zoom": { "7-9": 120, "10-11": 40, "12-14": 12 },
  "simplify_m_by_zoom": { "7-9": 80, "10-11": 25, "12-14": 8 },
  "drop_probability_by_zoom": { "7-9": 0.15, "10-11": 0.08, "12-14": 0.03 },
  "level_noise_m_by_zoom": { "7-9": 8, "10-11": 4, "12-14": 2 }
}
```

If `enabled=false`, the tool still reads settings but does not apply them.

### 8) Per-tile bathymetry payload options

This tool supports two broad strategies:

#### Strategy A — Embed bathy payload inside the vector tile (MVT)
Use the `bathy_data` block with `output="embedded_mvt"` and choose a mode.

```json
"bathy_data": {
  "enabled": true,
  "output": "embedded_mvt",
  "mode": "mvt_grid_i16",
  "layer_name": "bathy_data",
  "grid_size": 32,
  "quantization_m": 0.1
}
```

Supported embedded modes:
- `none` (stores nothing)
- `mvt_grid_i16` (zlib+base64 int16 grid)
- `mvt_grid_f32` (zlib+base64 float32 grid)
- `mvt_png` (zlib+base64 PNG bytes with vmin/vmax)

#### Strategy B — Write bathy payload to a separate **raster MBTiles** (PNG tiles)
This produces a second MBTiles containing PNG tiles aligned to the same z/x/y scheme.

```json
"bathy_data": {
  "enabled": true,
  "output": "raster_mbtiles",
  "raster_mbtiles": {
    "output_mbtiles": "./out/bathy_raster.mbtiles",
    "encoding": "png8",
    "grid_size": 256,
    "value_scaling": {
      "mode": "linear",
      "vmin": -6000,
      "vmax": 0
    }
  }
}
```

Raster MBTiles details:
- Each tile is a PNG image (typically 256×256)
- Values are encoded in pixel intensity (0..255) for `png8`
- The MBTiles `metadata` includes:
  - encoding
  - vmin/vmax (or other scaling parameters)
  - nodata conventions

Recommended uses:
- WebGL hillshading / local derivatives on the client
- Avoiding huge vector-tile properties when you want dense fields

Tradeoffs:
- Two tile sources to serve (vector contours + raster field)
- Raster is heavier than a small grid payload, lighter than full-resolution grids embedded in MVT

---

## Full config example (includes new options)

```json
{
  "dtm_year": 2024,
  "csw_url": "https://emodnet.ec.europa.eu/geonetwork/emodnet/eng/csw",

  "output_mbtiles": "./out/bathy_contours.mbtiles",
  "cache_dir": "./cache",

  "bbox_lonlat": [13.90, 40.55, 14.45, 40.95],
  "zoom_min": 7,
  "zoom_max": 14,

  "tile_pixels": 512,
  "extent": 4096,
  "layer_name": "bathy_contours",

  "contour_intervals_m": {
    "7-9": 50,
    "10-11": 20,
    "12-14": 10
  },

  "cache_intermediate": true,

  "coastline_mask": {
    "enabled": true,
    "source": "osm_land_polygons",
    "osm_land_polygons_url": "https://osmdata.openstreetmap.de/data/land-polygons-split-4326.zip"
  },

  "land_mask_cache": {
    "enabled": true,
    "all_touched": false
  },

  "bathy_data": {
    "enabled": true,
    "output": "embedded_mvt",
    "mode": "mvt_grid_i16",
    "layer_name": "bathy_data",
    "grid_size": 32,
    "quantization_m": 0.1,

    "raster_mbtiles": {
      "output_mbtiles": "./out/bathy_raster.mbtiles",
      "encoding": "png8",
      "grid_size": 256,
      "value_scaling": { "mode": "linear", "vmin": -6000, "vmax": 0 }
    }
  },

  "obfuscation": {
    "enabled": true,
    "seed": 20260201,
    "jitter_m_by_zoom": { "7-9": 120, "10-11": 40, "12-14": 12 },
    "simplify_m_by_zoom": { "7-9": 80, "10-11": 25, "12-14": 8 },
    "drop_probability_by_zoom": { "7-9": 0.15, "10-11": 0.08, "12-14": 0.03 },
    "level_noise_m_by_zoom": { "7-9": 8, "10-11": 4, "12-14": 2 }
  },

  "mbtiles_metadata": {
    "name": "Bathymetry contours (EMODnet-derived)",
    "description": "Contour tiles generated from EMODnet Bathymetry DTM tiles; land-masked and optionally obfuscated.",
    "attribution": "Derived from EMODnet Bathymetry DTM (CC BY 4.0)."
  }
}
```

---

## Client-side decoding notes

### Embedded grid (`mvt_grid_i16` / `mvt_grid_f32`)
- base64 decode `grid_b64`
- zlib decompress
- interpret bytes:
  - i16: int16 little-endian; depth_m = int16 * q_m (nodata = -32768)
  - f32: float32 little-endian; nodata = NaN
- reshape into (grid_h, grid_w)

### Raster MBTiles
- Fetch PNG tile (z/x/y) from the raster MBTiles endpoint
- Decode pixel intensity
- Map to depth using `vmin/vmax` from metadata (linear mode)

---

## Performance tips

- Keep `cache_intermediate=true`
- Keep `land_mask_cache.enabled=true` (big speed win)
- If bbox is small & zoom_max high, increase `tile_pixels` for smoother contours
- If bandwidth matters, prefer `bathy_data.output="embedded_mvt"` with small `grid_size` (e.g., 16 or 32)
- If GPU processing is the goal, prefer `bathy_data.output="raster_mbtiles"` and serve raster tiles separately

---

## Troubleshooting

### No tiles found
- Try another `dtm_year`
- Confirm bbox intersects EMODnet DTM coverage
- CSW endpoint may be temporarily down

### Fiona/Rasterio install issues
- Fiona/Rasterio depend on GDAL; use platform wheels where possible.
- On Linux, you may need system GDAL packages if wheels aren’t available.

### Empty output tiles
- Bbox is on land (everything masked)
- Bathymetry has nodata in that region (no contours possible)

---
