#!/usr/bin/env bash

# Assumes Deployer has ssh access to target instance and Provision scrippt has completed
scp docker-compose-prod.yaml ubuntu@global.infrastructureresilience.org:~/
scp -r envs/prod ubuntu@global.infrastructureresilience.org:/opt/infra-risk
ssh ubuntu@global.infrastructureresilience.org 'cd /opt/infra-risk && sudo chown -R gitaction:gitaction /opt/infra-risk && sudo su gitaction && cd /opt/infra-risk && docker-compose -f docker-compose-prod.yaml pull'
# Start DB
docker compose -f docker-compose-prod.yaml up -d db
# Setup DB schema using Backend image
docker compose -f docker-compose-prod.yaml run backend pipenv run python /code/create_schema.py
# Start other services
docker compose -f docker-compose-prod.yaml up -d backend raster-tilerserver vector-tileserver
# Start webserver
docker compose -f docker-compose-prod.yaml up -d web-server

# #
# # Deploy app
# #
# BASEDIR=$(dirname $0)

# pushd $BASEDIR/..

# echo "Running in ${pwd}"

# npm run build

# # built files for frontend
# rsync -rvz build/ ubuntu@jamaica.infrastructureresilience.org:/var/www/html

# # vector
# rsync -rvz tileserver/vector/data/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector/data
# rsync -rvz tileserver/vector/fonts/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector/fonts
# rsync -rvz tileserver/vector/config.json ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/vector

# # raster
# rsync -rvz tileserver/raster/data/ ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/raster/data
# rsync -rvz tileserver/raster/config.toml ubuntu@jamaica.infrastructureresilience.org:/var/www/tileserver/raster

# # backend
# rsync -rvz backend/backend/ ubuntu@jamaica.infrastructureresilience.org:/var/www/backend/backend
# rsync -rvz backend/setup.py ubuntu@jamaica.infrastructureresilience.org:/var/www/backend/
# pushd backend
#   pipenv lock -r > requirements.txt
# popd
# rsync -rvz backend/requirements.txt ubuntu@jamaica.infrastructureresilience.org:/var/www/backend/

# # update requirements
# ssh ubuntu@jamaica.infrastructureresilience.org 'cd /var/www/backend && sudo python3.10 -m pip install -r requirements.txt'
# # restart backend
# ssh ubuntu@jamaica.infrastructureresilience.org 'sudo service backend restart'
# ssh ubuntu@jamaica.infrastructureresilience.org 'sudo chown :www-data /var/www/backend/backend.sock'

# # restart backend
# ssh ubuntu@jamaica.infrastructureresilience.org 'sudo service terracotta restart'
# ssh ubuntu@jamaica.infrastructureresilience.org 'sudo chown :www-data /var/www/tileserver/raster/terracotta.sock'

# # restart tileserver
# ssh ubuntu@jamaica.infrastructureresilience.org 'sudo service tileservergl restart'

# popd
