#!/bin/env bash

RASTER_BASE_PATH=$1

terracotta ingest \
  "$RASTER_BASE_PATH/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__gcm_{gcm}.tif" \
  -o "$RASTER_BASE_PATH/terracotta.sqlite"
