# Deploy

The site can run on a single Linux virtual machine, with a separate database server.

The virtual machine runs several services, coordinated by docker-compose:

- Traefik proxy receives requests from the web, terminates SSL connections and
  handles basic authentication.
- Frontend React application, built using node and npm. In production this is
  stored and served as static files (HTML/JS/CSS).
- Vector tileserver, tileserver-gl-light, depends on node
- Raster tileserver, terracotta, depends on gunicorn and Python 3.10
- Backend Python application, depends on uvicorn, fastapi and Python 3.10

The frontend application source code is held in
[this repository](https://github.com/nismod/infra-risk-vis/) and this guide
assumes that it is built using node and npm locally on a development machine.
It would be possible to build directly on the server in a working directory.

To build and deploy the site:

- provision a server
- configure the server
- build and push docker images
- upload data and configuration
- load or restore data to the database
- pull and run the services

## AWS (optional)

This is optional, and only relevant if setting up on Amazon Web Services.

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

## On-premises (optional)

This is optional, and only relevant if setting up servers on-premises.

The service requires two virtual machines with similar specifications:

1. Application server
   - exposed to internal or public network to serve the app (ports 80 and 443)
   - running Ubuntu 20.04
   - minimum resources of ~20GB disk, ~1GB RAM, 2 cores
2. Database server
   - exposed to the application server
   - running Ubuntu 20.04
   - minimum resources of ~80GB disk, ~1GB RAM, 2 cores

## VM provisioning

`provision.sh` contains installation instructions for an Ubuntu 20.04 server to
install docker and docker-compose.

> This is relevant in either AWS or on-premises setup.

## Database provisioning

`provision-database-server.sh` contains installation instructions for an Ubuntu
20.04 server to install a PostgreSQL database server with the PostGIS extension.

> If running on AWS, this should be handled by terraform as an RDS database.

## Basic authentication

Create a password file for HTTP Basic Authentication if it doesn't already exist:

```bash
sudo touch /var/www/auth/.htpasswd
```

Add a user to the password file (will prompt for password):

```bash
sudo htpasswd -B /var/www/auth/.htpasswd new-username
```

## Database connection

The `PG*` variables for connection to database, to be replaced with actual details:

```conf
PGDATABASE=jamaica
PGUSER=docker
PGPASSWORD=docker
PGHOST=localhost
```

Testing database connection:

```bash
cd /var/www/
set -a
source ./envs/prod/.backend.env  # to use app connection details
# source ./envs/prod/.dbrestore.env  # to use admin connection details
set +a
# if needed:
# sudo apt install postgresql-client
psql
```

Restore a database dump:

```bash
pg_restore -cC -d postgres /path/to/backup.dump
```

## Deployment

Run `deploy.sh` to upload data and docker-compose config.

### Manage production

On remote, to first start the application:

```bash
cd /var/www
docker compose -f docker-compose.prod.yml up -d
```

### Publishing docker images

For pushing to the GitHub container registry, you will need to follow these
[instructions for authentication](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry).

Build and publish all images:

```bash
# Build images locally
docker compose -f docker-compose-prod.yaml build

# Push to GitHub container registry
docker push ghcr.io/nismod/jsrat-frontend:0.1
docker push ghcr.io/nismod/jsrat-backend:0.1
docker push ghcr.io/nismod/jsrat-vector-tileserver:0.1
docker push ghcr.io/nismod/jsrat-raster-tileserver:0.1
```

### Updating a service

Update a specific image, then build and push:

```bash
# Edit the image version in `docker-compose.prod.yaml`
# in this example it's on line 33:
#     image: ghcr.io/nismod/jsrat-frontend:0.1

# Build
docker compose -f docker-compose.prod.yaml build frontend

# Push
docker push ghcr.io/nismod/jsrat-frontend:0.1
```

Run `deploy.sh` to update the docker-compose config on the server.

On the remote server, pull the image, then restart the specific service:

```bash
# Pull image
docker pull ghcr.io/nismod/jsrat-frontend:0.1

# Restart service
docker compose up -d frontend
```
