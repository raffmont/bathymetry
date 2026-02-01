# obfuscate_bathy_tiles

A Python tool to generate obfuscated or exact bathymetric vector and raster tiles
from EMODnet Bathymetry DTM datasets.

## Features
- Automatic EMODnet DTM discovery & download
- On-the-fly raster mosaicking
- Per-zoom raster cache
- Coastline masking (OSM or EMODnet)
- Vector contours (MVT)
- Raster bathymetry (PNG MBTiles)
- Derived rasters (hillshade, slope)
- TileJSON generation
- Ready-to-use FastAPI demo server

See:
- [docs/maplibre_demo.md](maplibre_demo.md)
- [docs/openlayers_demo.md](openlayers_demo.md)
- [docs/fastapi_server.md](fastapi_server.md)
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
- `demo/maplibre.html`
- `demo/openlayers.html`

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

## Documentation index
- MapLibre: `docs/maplibre_demo.md`
- OpenLayers: `docs/openlayers_demo.md`
- Server: `docs/fastapi_server.md`
- Signal K usage: `docs/signalk_charts.md`
