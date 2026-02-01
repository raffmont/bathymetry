# OpenLayers demo

This demo renders **vector bathymetric contours (MVT)** and optional **raster overlays** (bathymetry PNG8 + hillshade PNG8)
produced by `obfuscate_bathy_tiles.py`.

## Files
- `demo/public/openlayers.html` – the demo page
- `demo/server.py` – a tiny FastAPI tile server that serves the exact endpoints expected by the demo
- `demo/config.json` – server configuration (tile locations + public directory)
- Tile outputs (generated):
  - `./data/output/bathy_contours.mbtiles` (vector, required)
  - `./data/output/bathy_raster.mbtiles` (optional)
  - `./data/output/bathy_hillshade.mbtiles` (optional)
  - `./data/output/bathy_slope.mbtiles` (optional)
  - `./data/output/tilejson/*.tilejson` (optional, for integration)

## 1) Generate tiles (once)
Example:
```bash
python obfuscate_bathy_tiles.py config.json
```

If you want **scientific (non-obfuscated)** contours:
```bash
python obfuscate_bathy_tiles.py config.json --no-obfuscation
```

## 2) Start the demo server
From the project root:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements.txt
python demo/server.py
```

The server defaults to:
- `http://127.0.0.1:8080`

Quick check:
- `http://127.0.0.1:8080/health`

### Custom MBTiles locations
Update `demo/config.json` if your outputs live elsewhere.

## 3) Open the demo page
Open `http://127.0.0.1:8080/openlayers.html` in your browser.

The page expects these endpoints:
- `/tiles/vector/{z}/{x}/{y}.pbf`
- `/tiles/bathy/{z}/{x}/{y}.png`
- `/tiles/hillshade/{z}/{x}/{y}.png`
- `/tiles/slope/{z}/{x}/{y}.png`

## Troubleshooting
- 404 on tiles: verify MBTiles paths and that you generated the zoom range you are viewing.
- No bathy/hillshade overlay: those MBTiles are optional; if missing you’ll see 404s (safe to ignore).
- CORS issues: `demo/server.py` enables CORS for `*`.
