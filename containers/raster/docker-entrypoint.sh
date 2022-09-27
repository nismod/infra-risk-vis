#!/bin/bash
set -e

if [ "$1" = 'serve' ]; then

  TC_DRIVER_PATH="/opt/terracotta/db/terracotta.sqllite"
  if [ ! -f "$TC_DRIVER_PATH" ]; then
    echo "Cannot find DB, aborting"
    exit -1
  fi
  TC_ALLOWED_ORIGINS_METADATA='["*"]'
  TC_ALLOWED_ORIGINS_TILES='["*"]'
  TC_PNG_COMPRESS_LEVEL=0
  TC_RESAMPLING_METHOD="nearest"
  TC_REPROJECTION_METHOD="nearest"
  gunicorn --workers 2 \
      --bind 0.0.0.0:5000 \
      --umask 007 \
      --access-logfile '-' \
      terracotta.server.app:app
  exit 0
fi

exec "$@"
