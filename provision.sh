#!/usr/bin/env bash

#
# Provision virtual machine
# - assuming OS is Ubuntu 18.04 LTS (Bionic)
#

# Install NGINX as proxy server
sudo apt-get update
sudo apt-get install nginx

# Set up SSL
sudo apt-get install software-properties-common
sudo add-apt-repository universe
sudo add-apt-repository ppa:certbot/certbot
sudo apt-get update
sudo apt-get install certbot python-certbot-nginx
sudo certbot certonly --nginx

# node version and platform
NODE_VERSION=v8.11.3
DISTRO=linux-x64

# download
wget -nc https://nodejs.org/dist/$NODE_VERSION/node-$NODE_VERSION-$DISTRO.tar.xz

# extract
sudo mkdir /usr/local/lib/node
sudo tar xf node-$NODE_VERSION-$DISTRO.tar.xz -C /usr/local/lib/node
rm node-$NODE_VERSION-$DISTRO.tar.xz

# setup
sudo mv /usr/local/lib/node/node-$NODE_VERSION-$DISTRO /usr/local/lib/node/node-$NODE_VERSION
sudo ln -s /usr/local/lib/node/node-$NODE_VERSION/bin/node /usr/bin/node

# install tileserver
npm i -g tileserver-gl-light
