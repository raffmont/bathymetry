# MapLibre demo

This demo renders **vector bathymetric contours (MVT)** and optional **raster overlays** (bathymetry PNG8 + hillshade PNG8)
produced by `obfuscate_bathy_tiles.py`.

## Files
- `demo/public/maplibre.html` – the demo page
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
Open `http://127.0.0.1:8080/maplibre.html` in your browser.

The page expects these endpoints:
- `/tiles/vector/{z}/{x}/{y}.pbf`
- `/tiles/bathy/{z}/{x}/{y}.png`
- `/tiles/hillshade/{z}/{x}/{y}.png`
- `/tiles/slope/{z}/{x}/{y}.png`

## Extending the MapLibre style

### A) Change contour styling
In `maplibre.html`, edit the `contours` layer:
- `"line-width"`
- `"line-opacity"`

If you want **depth-dependent styling**, you can use expressions on `depth_m`, e.g.
```js
"line-width": [
  "interpolate", ["linear"], ["abs", ["get", "depth_m"]],
  0, 0.3,
  50, 0.6,
  200, 1.2
]
```

### B) Colorize bathymetry raster
Current bathy is grayscale PNG8. You can:
- keep it as a background luminance
- switch to client-side color ramps with WebGL (custom layer)
- or output RGB tiles instead (customize generator)

### C) Add labels
You can add a symbol layer reading `depth_m` or `depth_label` as text.
Note: For dense contour lines, labeling needs decluttering logic and careful filters.

If you enable `"contour_labels"` in your generator config, each contour feature will include
`depth_label` formatted with your preferred string (see [docs/obfuscate_bathy_tiles.md](obfuscate_bathy_tiles.md)).

## Troubleshooting
- 404 on tiles: verify MBTiles paths and that you generated the zoom range you are viewing.
- No bathy/hillshade overlay: those MBTiles are optional; if missing you’ll see 404s (safe to ignore).
- CORS issues: `demo/server.py` enables CORS for `*`.
