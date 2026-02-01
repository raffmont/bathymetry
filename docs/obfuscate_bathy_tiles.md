# obfuscate_bathy_tiles

A Python tool to generate obfuscated or exact bathymetric vector and raster tiles
from EMODnet Bathymetry DTM datasets.

## Features
- Automatic EMODnet DTM discovery & download
- On-the-fly raster mosaicking (auto-reprojects mismatched CRS tiles to a common CRS)
- CRS inference from filenames that include EPSG codes (e.g., 3857 or 4326)
- Per-zoom raster cache
- GeoTIFF, ASCII grid, and NetCDF/HDF5 raster ingestion
- Multi-band rasters are ingested using the first band only
- Coastline masking (OSM or EMODnet)
- Vector contours (MVT)
- Raster bathymetry (PNG MBTiles)
- Derived rasters (hillshade, slope)
- TileJSON generation
- Ready-to-use FastAPI demo server

See:
- [docs/maplibre.md](maplibre.md)
- [docs/openlayers.md](openlayers.md)
- [docs/server.md](server.md)
- [docs/signalk_charts.md](signalk_charts.md)

---

## Demo & serving tiles locally (new)

### FastAPI demo server
A minimal server is provided in:
- `demo/server.py`

Install and run:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements.txt
python demo/server.py
```

Then open:
- `http://127.0.0.1:8080/maplibre.html`
- `http://127.0.0.1:8080/openlayers.html`

### Tile endpoints exposed
- `/tiles/vector/{z}/{x}/{y}.pbf`
- `/tiles/bathy/{z}/{x}/{y}.png`
- `/tiles/hillshade/{z}/{x}/{y}.png`
- `/tiles/slope/{z}/{x}/{y}.png`
- `/tilejson/{name}`

## New configuration keys (highlights)

### derived_rasters
Enable hillshade & slope MBTiles:
```json
"derived_rasters": {
  "enabled": true,
  "hillshade_mbtiles": "./data/output/bathy_hillshade.mbtiles",
  "slope_mbtiles": "./data/output/bathy_slope.mbtiles",
  "tile_size": 256,
  "hillshade": { "azimuth_deg": 315.0, "altitude_deg": 45.0 }
}
```

### tilejson
Write TileJSON files (recommended):
```json
"tilejson": {
  "enabled": true,
  "base_url": "http://localhost:8080/tiles",
  "out_dir": "./data/output/tilejson"
}
```

### CLI overrides
Run with config and override a few parameters without editing JSON:
```bash
python obfuscate_bathy_tiles.py config.json   --zoom-min 8 --zoom-max 13   --bbox 13.95 40.60 14.35 40.90   --no-obfuscation --verbose
```
Use `--verbose` to log CSW query details when troubleshooting tile discovery.

### contour_labels
Add a formatted label attribute to contour line features:
```json
"contour_labels": {
  "enabled": true,
  "format": "{depth_m} m"
}
```
The formatter accepts `{depth_m}` (rounded) and `{depth_m_raw}` (unrounded) placeholders.

### contour_smoothing
Smooth contours with Chaikin subdivision iterations per zoom:
```json
"contour_smoothing": {
  "iterations_by_zoom": {
    "7-9": 1,
    "10-11": 1,
    "12-14": 2
  }
}
```

## Documentation index
- MapLibre: `docs/maplibre.md`
- OpenLayers: `docs/openlayers.md`
- Server: `docs/server.md`
- Signal K usage: `docs/signalk_charts.md`
