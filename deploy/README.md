# Deploy

The site can run on a single Linux virtual machine.

The virtual machine runs several processes:
- Nginx reverse proxy and static file server, receives requests from the web,
  terminates SSL connections and handles basic authentication.
- Frontend React application, built using node and npm. On the server this is
  built into static files (HTML/JS/CSS).
- Vector tileserver, tileserver-gl-light, depends on nodejs
- Raster tileserver, terracotta, depends on gunicorn and Python

The frontend application source code is held in
[this repository](https://github.com/nismod/infra-risk-vis/) and this guide
assumes that it is built using node and npm locally on a development machine.
It would be possible to build directly on the server in a working directory.

To build and deploy the site:

- provision a server
- configure the server
- build and upload frontend, data and config

Server provision (and related DNS/access configuration) for AWS can be run using
[terraform](https://www.terraform.io/).

The scripts and configuration referenced below could be adapted to set up a
virtual machine or server in other environments, the AWS-specific elements are
all used to manage DNS and access.

Install the
[AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html)
then run:

```bash
aws configure   # one-off to set your AWS credentials
```

Install [terraform](https://www.terraform.io/) then run:

```bash
terraform init  # one-off to fetch provider from terraform registry
terraform plan  # to see what actions will be taken in detail
terraform apply # rerun after any change to main.tf
```

`provision.sh` contains installation instructions for an Ubuntu 20.04 server to
install NGINX, setup SSL using CertBot, install node and tileserver-gl-light

Create a password file for HTTP Basic Authentication:

```bash
sudo touch /etc/nginx/.htpasswd
```

Add a user to the password file (will prompt for password):

```bash
sudo htpasswd -B /etc/nginx/.htpasswd new-username
```

Configure terracotta to run under gunicorn:

```bash
sudo chown :www-data /var/www/tileserver/raster/terracotta.sock
sudo -u www-data curl --unix-socket /var/www/tileserver/raster/terracotta.sock http
```

Transfer tileserver raster data to the server, then ingest to terracotta:

```bash
terracotta ingest \
    "/var/www/tileserver/raster/data/{type}__rp_{rp}__rcp_{rcp}__epoch_{epoch}__conf_{confidence}.tif" \
    -o "/var/www/tileserver/raster/terracotta.sqlite"
```

`config/` directory contains:

- nginx config to serve frontend assets directly and proxy tile requests to the
  tileserver
- systemd service config to run the vector and raster tileservers as services

`deploy.sh` builds the frontend for deployment, uploads the build directory,
data and tileserver config to a server, and restarts the tileservers. It assumes
that whoever runs the script has ssh/public key access to the server.

## Deploy update

Usually this involves running `deploy.sh` and that will be sufficient.

If you need to restart the raster tileserver:

```bash
# restart the service
sudo service restart terracotta
# check nginx can access the socket
sudo chown :www-data /var/www/tileserver/raster/terracotta.sock
```
