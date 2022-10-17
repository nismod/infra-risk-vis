# Deploy

The site runs on a single Linux virtual machine using docker compose.

The stack consists of five services:
- Web server (nginx)
- Vector tileserver (tileserver-gl)
- Backend / API (bespoke Python app inc. raster tileserver via Terracotta Python API)
- Database (PostgreSQL with PostGIS)
- Tiles Database (MySQL)

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
install docker and docker compose.

`deploy.sh` uploads data to the server. It assumes that whoever runs the script
has ssh/public key access to the server.


#### Environment

The following env files are required:

##### .backend.env

```
PGHOST=
PGDATABASE=
PGUSER=
PGPASSWORD=

# Tiles API
LOG_LEVEL=INFO
RASTER_BASE_PATH=/data  # The mount underwich GeoTiffs for the tileserver can be found
MYSQL_URI=  # MySQL URI for tiles-db
API_TOKEN=  # Only required for mutating tiles metadata in the API

# Terracotta internal
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```


`docker compose up --build` will build the necessary images and bring up
containers with environment data provided from the .env file.

`docker ps` to see the containers' state.
