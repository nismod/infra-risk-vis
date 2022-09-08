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
rsync -rvz build/ ubuntu@east-africa.infrastructureresilience.org:/var/www/html

# vector
rsync -rvz tileserver/vector/data/ ubuntu@east-africa.infrastructureresilience.org:/var/www/tileserver/vector/data
rsync -rvz tileserver/vector/fonts/ ubuntu@east-africa.infrastructureresilience.org:/var/www/tileserver/vector/fonts
rsync -rvz tileserver/vector/config.json ubuntu@east-africa.infrastructureresilience.org:/var/www/tileserver/vector

# raster
rsync -rvz tileserver/raster/data/ ubuntu@east-africa.infrastructureresilience.org:/var/www/tileserver/raster/data
rsync -rvz tileserver/raster/config.toml ubuntu@east-africa.infrastructureresilience.org:/var/www/tileserver/raster

# backend
rsync -rvz backend/backend/ ubuntu@east-africa.infrastructureresilience.org:/var/www/backend/backend
rsync -rvz backend/setup.py ubuntu@east-africa.infrastructureresilience.org:/var/www/backend/
pushd backend
  pipenv lock -r > requirements.txt
popd
rsync -rvz backend/requirements.txt ubuntu@east-africa.infrastructureresilience.org:/var/www/backend/

# update requirements
ssh ubuntu@east-africa.infrastructureresilience.org 'cd /var/www/backend && sudo python3.10 -m pip install -r requirements.txt'
# restart backend
ssh ubuntu@east-africa.infrastructureresilience.org 'sudo service backend restart'
ssh ubuntu@east-africa.infrastructureresilience.org 'sudo chown :www-data /var/www/backend/backend.sock'

# restart raster tileserver
ssh ubuntu@east-africa.infrastructureresilience.org 'sudo service terracotta restart'
ssh ubuntu@east-africa.infrastructureresilience.org 'sudo chown :www-data /var/www/tileserver/raster/terracotta.sock'

# restart vector tileserver
ssh ubuntu@east-africa.infrastructureresilience.org 'sudo service tileservergl restart'

popd
