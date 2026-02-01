#!/usr/bin/env bash
set -e

echo "=== Bathymetry demo runner ==="

PYTHON_BIN=${PYTHON_BIN:-python3}
VENV_DIR=${VENV_DIR:-.venv}
SERVER_PORT=${PORT:-8080}

if ! command -v $PYTHON_BIN >/dev/null 2>&1; then
  echo "ERROR: python not found. Set PYTHON_BIN or install Python 3."
  exit 1
fi

if [ ! -d "$VENV_DIR" ]; then
  echo "[1/4] Creating virtual environment ($VENV_DIR)"
  $PYTHON_BIN -m venv "$VENV_DIR"
fi

echo "[2/4] Activating virtual environment"
# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"

echo "[3/4] Installing demo dependencies"
pip install --upgrade pip >/dev/null
pip install -r demo/requirements.txt

echo "[4/4] Starting FastAPI tile server"
echo "-----------------------------------------"
echo "Vector tiles:     http://127.0.0.1:${SERVER_PORT}/tiles/vector/{z}/{x}/{y}.pbf"
echo "Bathy raster:     http://127.0.0.1:${SERVER_PORT}/tiles/bathy/{z}/{x}/{y}.png"
echo "Hillshade raster: http://127.0.0.1:${SERVER_PORT}/tiles/hillshade/{z}/{x}/{y}.png"
echo "Slope raster:     http://127.0.0.1:${SERVER_PORT}/tiles/slope/{z}/{x}/{y}.png"
echo "Health:           http://127.0.0.1:${SERVER_PORT}/health"
echo ""
echo "Open these pages in your browser:"
echo "  http://127.0.0.1:${SERVER_PORT}/"
echo "  http://127.0.0.1:${SERVER_PORT}/maplibre.html"
echo "  http://127.0.0.1:${SERVER_PORT}/openlayers.html"
echo "-----------------------------------------"

export HOST=127.0.0.1
export PORT=$SERVER_PORT

exec $PYTHON_BIN demo/server.py
