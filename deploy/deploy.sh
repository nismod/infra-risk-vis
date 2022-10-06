#!/usr/bin/env bash

# Assumes Deployer has ssh access to target instance and Provision script has completed
scp docker-compose-prod.yaml ubuntu@global.infrastructureresilience.org:~/
scp -r envs/prod ubuntu@global.infrastructureresilience.org:~/
# ... Manually SSH in and move assets to /opt/infra-risk, then chown to gitaction
# sudo mv docker-compose-prod.yaml /opt/infra-risk/
# sudo mkdir /opt/infra-risk/envs
# sudo mv prod /opt/infra-risk/envs/
ssh ubuntu@global.infrastructureresilience.org 'cd /opt/infra-risk && sudo chown -R gitaction:gitaction /opt/infra-risk && sudo su gitaction && cd /opt/infra-risk && docker compose -f docker-compose-prod.yaml pull'
# Start DB
docker compose -f docker-compose-prod.yaml up -d db
# Setup DB schema using Backend image
docker compose -f docker-compose-prod.yaml run backend python /code/backend/create_schema.py
# Start other services
docker compose -f docker-compose-prod.yaml up -d backend raster-tilerserver vector-tileserver
# Start webserver
docker compose -f docker-compose-prod.yaml up -d web-server

