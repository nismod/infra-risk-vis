# Deploy

To build and deploy the site:

- provision a server
- configure the server
- build and upload frontend, data and config

Server provision (and related DNS/access configuration) for AWS can be run
using terraform.

Install the
[AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html)
then run:

  aws configure   # one-off to set your AWS credentials

Install [terraform](https://www.terraform.io/) then run:

  terraform init  # one-off to fetch provider from terraform registry
  terraform plan  # to see what actions will be taken in detail
  terraform apply # rerun after any change to main.tf

`provision.sh` contains installation instructions for an Ubuntu 20.04 server to
install NGINX, setup SSL using CertBot, install node and tileserver-gl-light

  sudo touch /etc/nginx/.htpasswd
  sudo htpasswd -B /etc/nginx/.htpasswd username

Getting terracotta running under gunicorn (TODO look at setting up gunicorn socket as service)

    sudo chown :www-data /var/www/tileserver/raster/terracotta.sock
    sudo -u www-data curl --unix-socket /var/www/tileserver/raster/terracotta.sock http

    terracotta ingest \
      "/var/www/tileserver/raster/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif" \
      -o "/var/www/tileserver/raster/terracotta.sqlite"


`config/` directory contains:

- nginx config to serve frontend assets directly and proxy tile requests to the
  tileserver
- systemd service config to run the vector and raster tileservers as services

`deploy.sh` builds the frontend for deployment, uploads the build directory,
data and tileserver config to a server, and restarts the tileservers. It assumes
that whoever runs the script has ssh/public key access.
