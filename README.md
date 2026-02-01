# Bathymetry

This project generates **bathymetric contour vector tiles** and optional **raster products**
(bathymetry field PNG8, hillshade PNG8, slope PNG8) from **EMODnet bathymetry DTM tiles**,
with optional **contour obfuscation** and **coastline-based land masking**.

The output is suitable for:
- MapLibre / OpenLayers web maps (MVT + PNG)
- local MBTiles workflows
- integration as a chart/map resource in Signal K ecosystems

## What you generate
- `data/output/bathy_contours.mbtiles` – vector tiles (MVT/PBF)
- `data/output/bathy_raster.mbtiles` – bathymetry raster tiles (PNG8) *(optional)*
- `data/output/bathy_hillshade.mbtiles` – hillshade tiles (PNG8) *(optional)*
- `data/output/bathy_slope.mbtiles` – slope tiles (PNG8) *(optional)*
- `data/output/tilejson/*.tilejson` – TileJSON descriptors *(optional but recommended)*
- `data/cache/` – cached intermediate rasters and per-zoom landmask caches

## Quick start
1) Generate tiles:
```bash
python obfuscate_bathy_tiles.py config.json
```
Add `--verbose` to include CSW discovery details in the logs.

2) Serve tiles locally and run demos:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements.txt
python demo/server.py
```
Then open:
- `demo/maplibre.html`
- `demo/openlayers.html`

## Documentation
- Generator documentation: [docs/obfuscate_bathy_tiles.md](docs/obfuscate_bathy_tiles.md)
- MapLibre demo: [docs/maplibre_demo.md](docs/maplibre_demo.md)
- OpenLayers demo: [docs/openlayers_demo.md](docs/openlayers_demo.md)
- FastAPI server: [docs/fastapi_server.md](docs/fastapi_server.md)
- Signal K integration: [docs/signalk_charts.md](docs/signalk_charts.md)

## Notes
- For public publishing, consider enabling contour obfuscation in the config (`obfuscation.enabled=true`) and tuning parameters by zoom.
- Raster outputs are optional; you can disable them if you only want contours.
