#!/usr/bin/env bash
set -e
set -x

#
# Deploy app data
# - assumes that SSH config is set up to connect to "jsrat1" host
#
BASEDIR=$(dirname $0)

pushd $BASEDIR/..

echo "Running in ${pwd}"

# vector data
rsync -ravz tileserver/vector/data/ jsrat1:/var/www/tileserver/vector/data
rsync -ravz tileserver/vector/fonts/ jsrat1:/var/www/tileserver/vector/fonts
rsync -ravz tileserver/vector/config.json jsrat1:/var/www/tileserver/vector

# raster data
rsync -ravz tileserver/raster/data/ jsrat1:/var/www/tileserver/raster/data
rsync -ravz tileserver/raster/config.toml jsrat1:/var/www/tileserver/raster

# docker compose configuration
rsync -avz docker-compose.prod.yml jsrat1:/var/www/
rsync -avz envs/prod jsrat1:/var/www/envs/

popd
