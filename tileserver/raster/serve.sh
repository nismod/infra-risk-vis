#!/bin/env bash

RASTER_BASE_PATH=$1

terracotta -c "$RASTER_BASE_PATH/config.toml" serve -d "$RASTER_BASE_PATH/terracotta.sqlite"  --port 5000