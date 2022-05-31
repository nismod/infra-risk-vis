#!/bin/env

set -e

# working directory for container
# this should still be set from the Dockerfile, but we'll cd here to be sure
cd $RASTER_BASE_PATH
RASTER_BASE_PATH=/opt/terracotta

# path to database
TC_DRIVER_PATH=$RASTER_BASE_PATH/db/terracotta.sqlite

# if there's no terracotta database, try and build one from the data files
if [ ! -f "$TC_DRIVER_PATH" ]; then
    echo "Cannot find DB, will try to ingest..."
    terracotta ingest \
      "$RASTER_BASE_PATH/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif" \
      -o "$TC_DRIVER_PATH"
fi

# development server
#terracotta -c "$RASTER_BASE_PATH/config.toml" serve -d "$TC_DRIVER_PATH" --port 5000 --allow-all-ips

# production deployment
# TC_ prepended env vars are terracotta config
# ... not sure how to cleanly provide these as a file

# for gunicorn settings, see
# https://docs.gunicorn.org/en/latest/settings.html#settings for gunicorn flags

TC_ALLOWED_ORIGINS_METADATA='["*"]' \
TC_ALLOWED_ORIGINS_TILES='["*"]' \
TC_PNG_COMPRESS_LEVEL=0 \
TC_RESAMPLING_METHOD="nearest" \
TC_REPROJECTION_METHOD="nearest" \
TC_DRIVER_PATH=$TC_DRIVER_PATH gunicorn \
    --workers 4 \
    --bind 0.0.0.0:5000 \
    --umask 007 \
    --access-logfile '-' \
    terracotta.server.app:app
