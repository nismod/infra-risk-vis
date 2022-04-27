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
sudo service terracotta restart
# check nginx can access the socket
sudo chown :www-data /var/www/tileserver/raster/terracotta.sock
```


## Deploy backend

Initial setup

```bash
# Python 3.10 on Ubuntu 20.04 needs to use PPA
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.10  python3.10-distutils
curl https://bootstrap.pypa.io/get-pip.py | sudo python3.10
# Install/upgrade Terracotta app requirements
sudo python3.10 -m pip install --upgrade cython
sudo python3.10 -m pip install --upgrade numpy
sudo python3.10 -m pip install --upgrade --no-binary rasterio rasterio==1.3a4
sudo python3.10 -m pip install --upgrade gunicorn terracotta[recommended]
sudo python3.10 -m pip install testresources

# Create backend working directory
sudo mkdir /var/www/backend
sudo chown -R :ubuntu /var/www/backend/
sudo chmod -R 775 /var/www/backend/
```

Run `deploy.sh` locally to copy over the backend package, including Python
package `requirements.txt`

```bash
# Install Python requirements
cd /var/www/backend
sudo python3.10 -m pip install -r requirements.txt

# Create backend service
sudo touch /etc/systemd/system/backend.service
# edit or copy from ./config/etc/systemd/system/backend.service
sudo systemctl start backend.service
# Check status
systemctl status backend
journalctl -u backend
```

Copy secret `PG*` variables for connection to database in EnvironmentFile as
configured in backend.service, something like:

```conf
PGDATABASE=jamaica
PGUSER=docker
PGPASSWORD=docker
PGHOST=localhost
```

To restart:

```bash
# restart the service
sudo service restart backend
# check nginx can access the socket
sudo chown :www-data /var/www/backend/backend.sock
```

Testing database connection

cd /var/www/backend/
set -a
source .env.prod
set +a
sudo apt install postgresql-client
psql
