# Demo guide

The demo bundles a small FastAPI server plus two frontend viewers (MapLibre and OpenLayers)
so you can preview vector contours and raster overlays produced by this project.

## What the demo needs
The demo expects MBTiles in `data/output/` and TileJSON in `data/output/tilejson/`. You can
produce exactly what the demo needs with:

```bash
python obfuscate_bathy_tiles.py examples/config-demo.json
```

The relevant outputs are:
- `data/output/bathy_contours.mbtiles` (vector contours)
- `data/output/bathy_raster.mbtiles` (bathymetry raster)
- `data/output/bathy_hillshade.mbtiles` (hillshade raster)
- `data/output/bathy_slope.mbtiles` (slope raster)
- `data/output/tilejson/*.tilejson`

## Start the demo server
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements.txt
python demo/server.py
```

The server reads `demo/config.json`, which points to the MBTiles above and configures the
public directory for the frontend.

## Open the viewers
- `http://127.0.0.1:8080/`
- `http://127.0.0.1:8080/maplibre.html`
- `http://127.0.0.1:8080/openlayers.html`

## Learn more
- Server details: [docs/server.md](server.md)
- MapLibre viewer: [docs/maplibre.md](maplibre.md)
- OpenLayers viewer: [docs/openlayers.md](openlayers.md)
