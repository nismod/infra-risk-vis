#!/bin/env bash

RASTER_BASE_PATH=$1

terracotta optimize-rasters \
  -o "$RASTER_BASE_PATH/data" \
  --overwrite \
  --reproject \
  --nproc -1 \
  --resampling-method nearest \
  "$RASTER_BASE_PATH/input/*.tif"

terracotta ingest \
  "$RASTER_BASE_PATH/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif" \
  -o "$RASTER_BASE_PATH/terracotta.sqlite"
