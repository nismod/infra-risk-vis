#!/bin/env bash

set -e
set -x
set -o errexit

# TODO validate arguments

INPUT_FILE=$1
OUTPUT_FILE=$2
OUTPUT_LAYER_NAME=$3
INPUT_LAYER_NAME=${4:-""}

rm -f tmp.json

ogr2ogr \
  -t_srs epsg:4326 \
  -f GeoJSONSeq \
  tmp.json \
  "$INPUT_FILE" \
  "$INPUT_LAYER_NAME"

tippecanoe \
  --output="$OUTPUT_FILE" \
  --layer="$OUTPUT_LAYER_NAME" \
  --force \
  ./tmp.json

rm tmp.json
