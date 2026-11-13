#!/usr/bin/env bash
# Download the static inputs into data_raw/: Natural Earth countries and cities, GHSL 2025 population grid.
# The road network itself comes from extract.py.
set -euo pipefail
RAW="$(cd "$(dirname "$0")/.." && pwd)/data_raw"
mkdir -p "$RAW"
cd "$RAW"

NE=https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson
for f in ne_10m_admin_0_countries ne_10m_populated_places_simple; do
  [ -f "$f.geojson" ] || curl -fL -o "$f.geojson" "$NE/$f.geojson"
done

GHSL=https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2025_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip
TIF=GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif
if [ ! -f "$TIF" ]; then
  curl -fL -o ghsl.zip "$GHSL"
  unzip -o ghsl.zip "$TIF"
  rm ghsl.zip
fi
ls -lh "$RAW"
