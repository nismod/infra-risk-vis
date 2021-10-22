#!/usr/bin/env bash

#
# Provision virtual machine
# - assuming OS is Ubuntu 20.04 LTS
#

# Install apt packages (NGINX, build requirements, GDAL)
sudo apt-get update
sudo apt-get install -y \
  nginx \
  software-properties-common \
  build-essential \
  gdal-bin \
  libgdal-dev \
  python3-pip \
  apache2-utils

# Set up SSL
sudo snap install core
sudo snap refresh core
sudo snap install --classic certbot
sudo ln -s /snap/bin/certbot /usr/bin/certbot

# TODO - this is interactive, prompts for email, agreement, domain
sudo certbot certonly --nginx

# Install node
NODE_VERSION=v14.18.1
DISTRO=linux-x64

wget -nc https://nodejs.org/dist/$NODE_VERSION/node-$NODE_VERSION-$DISTRO.tar.xz

sudo mkdir /usr/local/lib/node
sudo tar xf node-$NODE_VERSION-$DISTRO.tar.xz -C /usr/local/lib/node
rm node-$NODE_VERSION-$DISTRO.tar.xz

sudo mv /usr/local/lib/node/node-$NODE_VERSION-$DISTRO /usr/local/lib/node/node-$NODE_VERSION
sudo ln -s /usr/local/lib/node/node-$NODE_VERSION/bin/node /usr/bin/node
sudo ln -s /usr/local/lib/node/node-$NODE_VERSION/bin/npm /usr/bin/npm
sudo ln -s /usr/local/lib/node/node-$NODE_VERSION/bin/npx /usr/bin/npx

sudo chown :ubuntu /usr/local/lib/node/node-$NODE_VERSION/lib/node_modules/
sudo chmod 775 /usr/local/lib/node/node-$NODE_VERSION/lib/node_modules/
sudo chown :ubuntu /usr/local/lib/node/node-$NODE_VERSION/bin/
sudo chmod 775 /usr/local/lib/node/node-$NODE_VERSION/bin/

npm config set python /usr/bin/python3

# Install vector tileserver
npm i -g tileserver-gl-light

# Install raster tileserver
sudo pip install cython # must be available first to build dependencies
sudo pip install gunicorn terracotta[recommended]

sudo mkdir -p /var/www/tileserver/raster/data
sudo mkdir -p /var/www/tileserver/vector/data
sudo chown -R :ubuntu /var/www/
sudo chmod 775 /var/www/html/
sudo chmod 664 /var/www/html/*.html
sudo chmod -R  775 /var/www/tileserver/
