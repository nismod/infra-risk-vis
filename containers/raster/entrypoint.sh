#!/bin/env

set -e

RASTER_BASE_PATH=/opt/terracotta  # working directory for container
DB_PATH=$RASTER_BASE_PATH/db/terracotta.sqlite

if [ ! -f "$DB_PATH" ]; then
    echo "Cannot find DB, will try to ingest..."
    terracotta ingest \
      "$RASTER_BASE_PATH/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif" \
      -o "$DB_PATH"
fi

terracotta -c "$RASTER_BASE_PATH/config.toml" serve -d "$DB_PATH" --port 5000
