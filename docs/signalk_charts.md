# Signal K integration: using generated tiles as a chart/map resource

This document describes how to expose the generated bathymetry products (vector contours + raster bathy/hillshade)
to a Signal K ecosystem as map resources.

## What you have
Generated outputs (typical):
- Vector contours: `bathy_contours.mbtiles` (MVT)
- Optional raster bathy: `bathy_raster.mbtiles` (PNG8)
- Optional hillshade: `bathy_hillshade.mbtiles` (PNG8)
- Optional slope: `bathy_slope.mbtiles` (PNG8)

## Serving tiles
Signal K clients typically need HTTP endpoints with XYZ templates.

Using the provided FastAPI server, you already have:
- `/tiles/vector/{z}/{x}/{y}.pbf`
- `/tiles/bathy/{z}/{x}/{y}.png`
- `/tiles/hillshade/{z}/{x}/{y}.png`
- `/tiles/slope/{z}/{x}/{y}.png`

Run it on the same machine as Signal K (or reachable from it).

## Registering as a map resource
Exact steps depend on the Signal K server UI/plugin you use.
Common approaches:

### A) Add as a custom chart layer in your client
Many Signal K chart clients (including those embedding MapLibre/OpenLayers) allow adding custom XYZ sources.
Add a vector layer for contours and (optionally) raster overlays.

Example templates:
- Vector: `http://<host>:8080/tiles/vector/{z}/{x}/{y}.pbf`
- Raster: `http://<host>:8080/tiles/bathy/{z}/{x}/{y}.png`
- Hillshade: `http://<host>:8080/tiles/hillshade/{z}/{x}/{y}.png`

### B) Use TileJSON (recommended)
If your client supports TileJSON, serve:
- `http://<host>:8080/tilejson/bathy_contours.tilejson`
- `http://<host>:8080/tilejson/bathy_raster.tilejson`
- `http://<host>:8080/tilejson/hillshade.tilejson`
- `http://<host>:8080/tilejson/slope.tilejson`

TileJSON provides bounds, zooms, and the URL templates.

## Styling in nautical contexts
Suggested approach:
- Keep bathymetry raster faint (opacity ~0.3–0.6)
- Use hillshade to give seabed shape perception (opacity ~0.2–0.5)
- Use contours for legibility and depth awareness
- If you need safety/obfuscation, enable `obfuscation.enabled=true` and tune by zoom

## Client-side depth decoding (PNG8)
If you use raster bathy PNG8:
- 0 pixel value = nodata (land or missing)
- depths are linearly mapped between `vmin..vmax` (config: `bathy_data.raster_mbtiles.value_scaling`)
- decode: `depth = vmin + (pixel/255)*(vmax-vmin)`

## Security / obfuscation notes
For public layers, consider:
- enabling obfuscation for contours
- using coarser contour intervals at high zoom
- reducing raster tile resolution or disabling raster output

## Extending
You can add additional derived products:
- curvature / rugosity
- aspect
- gradient magnitude
and serve them as `/tiles/<product>/{z}/{x}/{y}.png` similarly.
