#!/bin/env bash

set -e
set -x
set -o errexit

# TODO validate arguments

INPUT_FILE=$1
OUTPUT_FILE=$2
OUTPUT_LAYER_NAME=$3
SPATIAL_TYPE=$4
INPUT_LAYER_NAME=${5:-""}
INPUT_FILTER=${6:-"1=1"}

if [ $SPATIAL_TYPE == 'line' ]; then
  TIPPECANOE_OPTIONS="--drop-smallest-as-needed --minimum-zoom=3 --maximum-zoom=15"
elif [ $SPATIAL_TYPE == 'polygon' ]; then
  TIPPECANOE_OPTIONS="--drop-smallest-as-needed --minimum-zoom=3 --maximum-zoom=15"
elif [ $SPATIAL_TYPE == 'point' ]; then
  TIPPECANOE_OPTIONS="-zg"
else
  echo "Unknown spatial type"
  exit 1
fi

rm -f tmp.json

ogr2ogr \
  -t_srs epsg:4326 \
  -f GeoJSONSeq \
  -where "$INPUT_FILTER" \
  tmp.json \
  "$INPUT_FILE" \
  "$INPUT_LAYER_NAME"

tippecanoe \
  --generate-ids \
  --read-parallel \
  --output="$OUTPUT_FILE" \
  --layer="$OUTPUT_LAYER_NAME" \
  $TIPPECANOE_OPTIONS \
  --force \
  ./tmp.json

rm tmp.json
