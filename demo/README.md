# Demo (MapLibre + OpenLayers)

These HTML files are static demos that assume you are serving your MBTiles through a tile server.

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
Serve them as static files if you like.

## Serving tiles
Use any MBTiles-capable server (custom FastAPI/Flask, Tegola, Martin, TileServer-GL, etc.).
