# FastAPI tile server

`demo/server.py` is a tiny **read-only** tile server that reads MBTiles (SQLite) and serves tiles
with the exact URL schema required by the MapLibre/OpenLayers demos.

## Endpoints
- `GET /health` – quick diagnostics
- `GET /tiles/vector/{z}/{x}/{y}.pbf` – vector tiles (MVT)
- `GET /tiles/bathy/{z}/{x}/{y}.png` – bathy raster tiles (PNG8)
- `GET /tiles/hillshade/{z}/{x}/{y}.png` – hillshade tiles (PNG8)
- `GET /tiles/slope/{z}/{x}/{y}.png` – slope tiles (PNG8)
- `GET /tilejson/{name}` – serves `*.tilejson` from the configured tilejson directory

## Run
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements-demo.txt
python demo/server.py
```

Defaults:
- host: `127.0.0.1`
- port: `8080`

Override with env vars:
```bash
HOST=0.0.0.0 PORT=8080 python demo/server.py
```

## Configure MBTiles paths
Set any of:
```bash
export BATHY_VECTOR_MBTILES=./out/bathy_contours.mbtiles
export BATHY_RASTER_MBTILES=./out/bathy_raster.mbtiles
export BATHY_HILLSHADE_MBTILES=./out/bathy_hillshade.mbtiles
export BATHY_SLOPE_MBTILES=./out/bathy_slope.mbtiles
export TILEJSON_DIR=./out/tilejson
```

## Implementation notes
- MBTiles store tiles using **TMS Y indexing**.
- Web clients request tiles using **XYZ Y indexing**.
- The server converts `y_xyz` to `y_tms` with: `(2^z - 1) - y_xyz`.
- Vector tiles are gzip-compressed on the fly and served with `Content-Encoding: gzip`.

## Extending for production
Consider adding:
- stronger caching headers (ETag), CDN, and HTTP/2
- request logging and metrics
- range support and/or streaming
- precompressed vector tiles (store gzipped blobs)
- authentication/authorization for restricted datasets
