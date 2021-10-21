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
rsync -rvz build/ ubuntu@seasia.infrastructureresilience.org:/var/www/html

# data and config for tileserver
rsync -rvz data/ ubuntu@seasia.infrastructureresilience.org:/var/www/tileserver/data
rsync -rvz styles/ ubuntu@seasia.infrastructureresilience.org:/var/www/tileserver/styles
rsync -rvz fonts/ ubuntu@seasia.infrastructureresilience.org:/var/www/tileserver/fonts
rsync -rvz config.json ubuntu@seasia.infrastructureresilience.org:/var/www/tileserver

# restart tileserver
ssh ubuntu@seasia.infrastructureresilience.org 'sudo service tileserver restart'

popd
