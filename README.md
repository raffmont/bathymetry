# Bathymetry

This project generates **bathymetric contour vector tiles** and optional **raster products**
(bathymetry field PNG8, hillshade PNG8, slope PNG8) from **EMODnet bathymetry DTM tiles**,
with optional **contour obfuscation** and **coastline-based land masking**.

The output is suitable for:
- MapLibre / OpenLayers web maps (MVT + PNG)
- local MBTiles workflows
- integration as a chart/map resource in Signal K ecosystems

## What you generate
- `out/bathy_contours.mbtiles` – vector tiles (MVT/PBF)
- `out/bathy_raster.mbtiles` – bathymetry raster tiles (PNG8) *(optional)*
- `out/bathy_hillshade.mbtiles` – hillshade tiles (PNG8) *(optional)*
- `out/bathy_slope.mbtiles` – slope tiles (PNG8) *(optional)*
- `out/tilejson/*.tilejson` – TileJSON descriptors *(optional but recommended)*
- `cache/` – cached intermediate rasters and per-zoom landmask caches

## Quick start
1) Generate tiles:
```bash
python obfuscate_bathy_tiles.py config.all_features.json
```

2) Serve tiles locally and run demos:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r demo/requirements-demo.txt
python demo/server.py
```
Then open:
- `demo/maplibre.html`
- `demo/openlayers.html`

## Documentation
- Generator documentation: `obfuscate_bathy_tiles.md`
- MapLibre demo: `docs/maplibre_demo.md`
- OpenLayers demo: `docs/openlayers_demo.md`
- FastAPI server: `docs/fastapi_server.md`
- Signal K integration: `docs/signalk_charts.md`

## Notes
- For public publishing, consider enabling contour obfuscation in the config (`obfuscation.enabled=true`) and tuning parameters by zoom.
- Raster outputs are optional; you can disable them if you only want contours.
