#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tiny FastAPI tile server for the Bathymetry demo."""  # Provide a short module description.

from __future__ import annotations  # Enable forward-referenced type annotations.

import gzip  # Compress vector tiles for efficient transfer.
import json  # Parse JSON configuration files.
import logging  # Emit structured logs instead of print statements.
import os  # Access environment variables for configuration overrides.
import sqlite3  # Read MBTiles (SQLite) databases.
from functools import lru_cache  # Reuse SQLite connections safely.
from pathlib import Path  # Handle filesystem paths in a cross-platform way.
from typing import Any, Dict, Optional  # Define type hints for clarity.

from fastapi import FastAPI, HTTPException, Response  # Build API endpoints and responses.
from fastapi.middleware.cors import CORSMiddleware  # Allow browser-based demos to fetch tiles.
from fastapi.staticfiles import StaticFiles  # Serve the demo HTML/CSS/JS assets.

LOGGER = logging.getLogger("demo.server")  # Create a module-specific logger.


def _configure_logging() -> None:
    """Configure basic logging for the demo server."""  # Explain the helper purpose.
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")  # Set default log format.


_configure_logging()  # Initialize logging for module-level startup tasks.


def _config_path() -> Path:
    """Resolve the JSON config path, optionally overridden by DEMO_CONFIG."""  # Document path resolution.
    return Path(os.environ.get("DEMO_CONFIG", "config.json")).expanduser().resolve()  # Load default config.


def _load_config(path: Path) -> Dict[str, Any]:
    """Load the JSON config file and return it as a dictionary."""  # Explain function return.
    if not path.exists():  # Ensure the config file exists before reading.
        raise FileNotFoundError(f"Config not found: {path}")  # Fail fast with a clear error.
    LOGGER.info("Loading config from %s", path)  # Log which config file is used.
    return json.loads(path.read_text(encoding="utf-8"))  # Parse the JSON text content.


def _resolve_path(value: str, base_dir: Path) -> Path:
    """Resolve a config path relative to the config directory."""  # Clarify resolution strategy.
    candidate = Path(value).expanduser()  # Expand user tildes before resolving the path.
    if candidate.is_absolute():  # Return immediately for already-absolute paths.
        return candidate.resolve()  # Resolve the absolute path for consistency.
    return (base_dir / candidate).resolve()  # Resolve relative paths against the config directory.


def _get_nested(config: Dict[str, Any], *keys: str, default: Any = None) -> Any:
    """Read nested configuration values with a default fallback."""  # Explain nested lookup helper.
    current: Any = config  # Start at the root config dictionary.
    for key in keys:  # Walk down the nested key path.
        if not isinstance(current, dict) or key not in current:  # Guard against missing keys.
            return default  # Return default when a key is missing.
        current = current[key]  # Move deeper into the config structure.
    return current  # Return the final nested value.


_config_path_value = _config_path()  # Resolve the config path at import time.
_config = _load_config(_config_path_value)  # Load the JSON configuration into memory.
_config_dir = _config_path_value.parent  # Anchor relative paths to the config file location.

VECTOR_MBTILES = _resolve_path(_get_nested(_config, "tiles", "vector_mbtiles", default="./data/output/bathy_contours.mbtiles"), _config_dir)  # Vector MBTiles path.
BATHY_MBTILES = _resolve_path(_get_nested(_config, "tiles", "bathy_mbtiles", default="./data/output/bathy_raster.mbtiles"), _config_dir)  # Bathy raster MBTiles path.
HS_MBTILES = _resolve_path(_get_nested(_config, "tiles", "hillshade_mbtiles", default="./data/output/bathy_hillshade.mbtiles"), _config_dir)  # Hillshade MBTiles path.
SLOPE_MBTILES = _resolve_path(_get_nested(_config, "tiles", "slope_mbtiles", default="./data/output/bathy_slope.mbtiles"), _config_dir)  # Slope MBTiles path.
TILEJSON_DIR = _resolve_path(_get_nested(_config, "tiles", "tilejson_dir", default="./data/output/tilejson"), _config_dir)  # TileJSON output directory.
PUBLIC_DIR = _resolve_path(_get_nested(_config, "public_dir", default="./demo/public"), _config_dir)  # Directory holding demo assets.

APP = FastAPI(title="Bathymetry MBTiles Server", version="0.2.0")  # Instantiate the FastAPI application.

APP.add_middleware(
    CORSMiddleware,  # Insert CORS middleware into the stack.
    allow_origins=["*"],  # Allow any origin for local demo convenience.
    allow_methods=["GET"],  # Limit CORS to GET requests.
    allow_headers=["*"],  # Permit all headers in demo requests.
)  # End CORS middleware configuration.


def _tms_y(z: int, y_xyz: int) -> int:
    """Convert XYZ Y coordinates to TMS Y coordinates."""  # Explain the conversion formula.
    return (2 ** z - 1) - y_xyz  # Compute the TMS row index for MBTiles.


@lru_cache(maxsize=8)  # Cache a handful of SQLite connections.
def _conn(path: str) -> sqlite3.Connection:
    """Open (or reuse) an SQLite connection for an MBTiles file."""  # Explain connection caching.
    con = sqlite3.connect(path, check_same_thread=False)  # Create a thread-safe SQLite connection.
    con.row_factory = sqlite3.Row  # Return rows as dict-like objects.
    return con  # Provide the cached connection.


def _fetch_tile(mbtiles: Path, z: int, x: int, y_xyz: int) -> Optional[bytes]:
    """Read a tile blob from an MBTiles file, or return None."""  # Describe the function behavior.
    if not mbtiles.exists():  # Skip if the MBTiles file is missing.
        return None  # Indicate absence to the caller.
    con = _conn(str(mbtiles))  # Get a cached SQLite connection.
    y = _tms_y(z, y_xyz)  # Convert XYZ Y to TMS Y for MBTiles lookup.
    row = con.execute(
        "SELECT tile_data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?",  # Query tile data.
        (z, x, y),  # Bind parameters for the tile coordinates.
    ).fetchone()  # Fetch the first matching tile row.
    return row["tile_data"] if row else None  # Return the tile blob or None.


def _pbf_response(data: bytes) -> Response:
    """Build a gzipped response for vector (PBF) tiles."""  # Note response type.
    gz = gzip.compress(data, compresslevel=6)  # Compress the vector tile payload.
    return Response(
        content=gz,  # Use gzipped data as the response body.
        media_type="application/x-protobuf",  # Use the PBF media type.
        headers={
            "Content-Encoding": "gzip",  # Tell clients the data is gzipped.
            "Cache-Control": "public, max-age=3600",  # Cache vector tiles for one hour.
        },
    )  # Return the configured response object.


def _png_response(data: bytes) -> Response:
    """Build a response for PNG tiles."""  # Explain the response intent.
    return Response(
        content=data,  # Use raw PNG data as the response body.
        media_type="image/png",  # Use the PNG media type.
        headers={"Cache-Control": "public, max-age=3600"},  # Cache raster tiles for one hour.
    )  # Return the configured response object.


@APP.get("/health")  # Expose a health endpoint for diagnostics.
def health() -> Dict[str, Any]:
    """Return diagnostics about configured tile sources."""  # Describe the health endpoint.
    return {
        "vector_mbtiles": str(VECTOR_MBTILES),  # Report the vector MBTiles path.
        "bathy_mbtiles": str(BATHY_MBTILES),  # Report the bathy MBTiles path.
        "hillshade_mbtiles": str(HS_MBTILES),  # Report the hillshade MBTiles path.
        "slope_mbtiles": str(SLOPE_MBTILES),  # Report the slope MBTiles path.
        "tilejson_dir": str(TILEJSON_DIR),  # Report the TileJSON directory.
        "public_dir": str(PUBLIC_DIR),  # Report the static assets directory.
        "vector_exists": VECTOR_MBTILES.exists(),  # Indicate if vector MBTiles exist.
        "bathy_exists": BATHY_MBTILES.exists(),  # Indicate if bathy MBTiles exist.
        "hillshade_exists": HS_MBTILES.exists(),  # Indicate if hillshade MBTiles exist.
        "slope_exists": SLOPE_MBTILES.exists(),  # Indicate if slope MBTiles exist.
        "tilejson_exists": TILEJSON_DIR.exists(),  # Indicate if TileJSON directory exists.
        "public_exists": PUBLIC_DIR.exists(),  # Indicate if public directory exists.
    }  # Return the diagnostics payload.


@APP.get("/tiles/vector/{z}/{x}/{y}.pbf")  # Register the vector tiles route.
def vector_tile(z: int, x: int, y: int) -> Response:
    """Serve vector contour tiles from the configured MBTiles file."""  # Describe the endpoint.
    blob = _fetch_tile(VECTOR_MBTILES, z, x, y)  # Attempt to read the vector tile.
    if blob is None:  # Check for missing data.
        raise HTTPException(status_code=404, detail="Vector tile not found (MBTiles missing or empty).")  # Signal 404.
    return _pbf_response(blob)  # Return a gzipped PBF response.


@APP.get("/tiles/bathy/{z}/{x}/{y}.png")  # Register the bathymetry raster route.
def bathy_tile(z: int, x: int, y: int) -> Response:
    """Serve bathymetry raster tiles from the configured MBTiles file."""  # Describe the endpoint.
    blob = _fetch_tile(BATHY_MBTILES, z, x, y)  # Attempt to read the bathy tile.
    if blob is None:  # Check for missing data.
        raise HTTPException(status_code=404, detail="Bathy raster tile not found (MBTiles missing or empty).")  # Signal 404.
    return _png_response(blob)  # Return a PNG response.


@APP.get("/tiles/hillshade/{z}/{x}/{y}.png")  # Register the hillshade raster route.
def hillshade_tile(z: int, x: int, y: int) -> Response:
    """Serve hillshade raster tiles from the configured MBTiles file."""  # Describe the endpoint.
    blob = _fetch_tile(HS_MBTILES, z, x, y)  # Attempt to read the hillshade tile.
    if blob is None:  # Check for missing data.
        raise HTTPException(status_code=404, detail="Hillshade tile not found (MBTiles missing or empty).")  # Signal 404.
    return _png_response(blob)  # Return a PNG response.


@APP.get("/tiles/slope/{z}/{x}/{y}.png")  # Register the slope raster route.
def slope_tile(z: int, x: int, y: int) -> Response:
    """Serve slope raster tiles from the configured MBTiles file."""  # Describe the endpoint.
    blob = _fetch_tile(SLOPE_MBTILES, z, x, y)  # Attempt to read the slope tile.
    if blob is None:  # Check for missing data.
        raise HTTPException(status_code=404, detail="Slope tile not found (MBTiles missing or empty).")  # Signal 404.
    return _png_response(blob)  # Return a PNG response.


@APP.get("/tilejson/{name}")  # Register the TileJSON route.
def tilejson(name: str) -> Response:
    """Serve TileJSON files from the configured directory."""  # Describe the endpoint.
    p = (TILEJSON_DIR / name).resolve()  # Build an absolute path to the requested TileJSON.
    if not str(p).startswith(str(TILEJSON_DIR)):  # Guard against path traversal.
        raise HTTPException(status_code=400, detail="Invalid name")  # Reject invalid names.
    if not p.exists():  # Ensure the TileJSON file exists.
        raise HTTPException(status_code=404, detail="TileJSON not found")  # Signal 404 if missing.
    return Response(content=p.read_bytes(), media_type="application/json", headers={"Cache-Control": "no-cache"})  # Serve JSON.


APP.mount("/", StaticFiles(directory=str(PUBLIC_DIR), html=True), name="public")  # Serve the demo UI assets.


if __name__ == "__main__":
    _configure_logging()  # Initialize logging when running directly.
    host = _get_nested(_config, "server", "host", default="127.0.0.1")  # Read the host from config.
    port = int(_get_nested(_config, "server", "port", default=8080))  # Read the port from config.
    reload_flag = bool(_get_nested(_config, "server", "reload", default=True))  # Read reload flag from config.
    LOGGER.info("Starting demo server on %s:%s", host, port)  # Log server startup details.
    import uvicorn  # Import uvicorn only when needed for running the server.
    uvicorn.run("server:APP", host=host, port=port, reload=reload_flag)  # Start the ASGI server.
