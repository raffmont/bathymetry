# Demo (MapLibre + OpenLayers)

These HTML files are static demos served by the FastAPI tile server.

## Expected endpoints

Vector (MVT):
- http://localhost:8080/tiles/vector/{z}/{x}/{y}.pbf

Raster bathymetry (PNG8):
- http://localhost:8080/tiles/bathy/{z}/{x}/{y}.png

Derived rasters (PNG8):
- http://localhost:8080/tiles/hillshade/{z}/{x}/{y}.png
- http://localhost:8080/tiles/slope/{z}/{x}/{y}.png

## TileJSON
The generator writes TileJSON files to `./data/output/tilejson/` (configurable via `tilejson.out_dir`).
Serve them from the demo server at `/tilejson/{name}`.

## Serving tiles
Use the FastAPI demo server:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements.txt
python demo/server.py
```

The server reads `demo/config.json` for MBTiles locations and serves the static demo UI from
`demo/public/`.
