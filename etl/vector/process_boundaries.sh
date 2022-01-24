#!/bin/env bash

set -e
set -x
set -o errexit

# TODO validate arguments

INPUT_FILE=$1
POLYGON_OUTPUT_FILE=$2
LABEL_OUTPUT_FILE=$3
OUTPUT_LAYER_NAME=$4
INPUT_LAYER_NAME=${5:-""}

rm -f polygons.json labels.json

ogr2ogr \
  -t_srs epsg:4326 \
  -f GeoJSONSeq \
  polygons.json \
  "$INPUT_FILE" \
  "$INPUT_LAYER_NAME"

cat polygons.json | geojson-polygon-labels > labels.json \
  --ndjson \
  --label=polylabel \
  --style=largest \
  --include-minzoom=9-16 \
  --include-area \
  --include-bbox

tippecanoe \
  --read-parallel \
  --output="$POLYGON_OUTPUT_FILE" \
  --force \
  polygons.json

tippecanoe \
  --read-parallel \
  --output="$LABEL_OUTPUT_FILE" \
  --force \
  labels.json

rm polygons.json labels.json
