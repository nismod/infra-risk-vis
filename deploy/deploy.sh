#!/usr/bin/env bash
set -e
set -x

#
# Deploy app
#
BASEDIR=$(dirname $0)

pushd $BASEDIR/..

echo "Running in ${pwd}"

npm run build

# built files for frontend
rsync -rvz build/ ubuntu@jamaica.infrastructureresilience.org:/var/www/html

# vector
rsync -rvz tileserver/vector/data/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector/data
rsync -rvz tileserver/vector/fonts/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector/fonts
rsync -rvz tileserver/vector/config.json ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector

# raster
rsync -rvz tileserver/raster/data/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/raster/data
rsync -rvz tileserver/raster/terracotta.sqlite ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/raster
rsync -rvz tileserver/raster/config.toml ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/raster

# restart tileserver
ssh ubuntu@jamaica.infrastructureresilience.org 'sudo service tileservergl restart'

popd
