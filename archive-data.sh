#!/bin/bash
set -e
set -x
DATESTR=$(date -I)
tar cJf "${DATESTR}_tileserver-vector-data.tar.xz" tileserver/vector/data/
tar cJf "${DATESTR}_tileserver-raster-data.tar.xz" tileserver/raster/data/
pg_dump -Fc $PGDATABASE > "${DATESTR}_${PGDATABASE}.dump"
