#!/usr/bin/env bash

#
# Provision virtual machine
# - assuming OS is Ubuntu 20.04 LTS
# - assuming this script is run as a user in groups sudo and jsrat_admin
#

# Install helper apt packages
sudo apt-get update
sudo apt-get install -y \
  ca-certificates curl gnupg \
  apache2-utils

# Install docker
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

# Set up current user with docker
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker

# Set up data directories
sudo mkdir -p /var/www/tileserver/raster/data
sudo mkdir -p /var/www/tileserver/vector/data
sudo chown -R :jsrat_admin /var/www/
sudo chmod -R  775 /var/www/tileserver/
