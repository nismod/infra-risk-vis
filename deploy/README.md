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

`provision.sh` contains installation instructions for an Ubuntu 22.04 server to
install docker and docker compose.

`deploy.sh` uploads data to the server. It assumes that whoever runs the script
has ssh/public key access to the server.

## Environment

The following env files are required:

### .backend.env

```
PYTHONPATH=/code/backend

PGHOST=
PGDATABASE=
PGUSER=
PGPASSWORD=

# Tiles API
LOG_LEVEL=INFO
RASTER_BASE_PATH=/data  # The mount under which GeoTIFFs for the tileserver can be found
MYSQL_URI=  # MySQL URI for tiles-db
API_TOKEN=  # Only required for mutating tiles metadata in the API
DOMAIN_TO_DB_MAP='{"land_cover":"land_cover", "fluvial": "aqueduct", "coastal": "aqueduct", "extreme_heat": "extreme_heat", "cyclone": "cyclone", "population": "jrc_pop", "earthquake": "gem_earthquake", "nature": "exposure_nature", "buildings":"buildings", "drought":"drought", "traveltime_to_healthcare":"traveltime_to_healthcare"}'


# Terracotta internal
TC_ALLOWED_ORIGINS_METADATA='["*"]'
TC_ALLOWED_ORIGINS_TILES='["*"]'
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"
```

### .irv-autopkg.env

```
# Setup Python Path
PYTHONPATH="/usr/src/app"
AUTOPKG_DEPLOYMENT_ENV=prod
AUTOPKG_STORAGE_BACKEND=awss3
AUTOPKG_LOCALFS_STORAGE_BACKEND_ROOT=/data/packages
AUTOPKG_LOCALFS_PROCESSING_BACKEND_ROOT=/data/processing
AUTOPKG_LOCALFS_STORAGE_BACKEND_ROOT_TEST=/usr/src/app/tests/data/packages
AUTOPKG_LOCALFS_PROCESSING_BACKEND_ROOT_TEST=/usr/src/app/tests/data/processing
AUTOPKG_LOG_LEVEL=INFO
AUTOPKG_POSTGRES_USER=
AUTOPKG_POSTGRES_PASSWORD=
AUTOPKG_POSTGRES_HOST=
AUTOPKG_POSTGRES_PORT=
AUTOPKG_POSTGRES_DB=
AUTOPKG_CELERY_BROKER=redis://redis:6379
AUTOPKG_CELERY_BACKEND=redis://redis:6379
AUTOPKG_CELERY_CONCURRENCY=4
AUTOPKG_TASK_LOCK_TIMEOUT=43200
AUTOPKG_TASK_EXPIRY_SECS=43200
AUTOPKG_REDIS_HOST=redis
#AUTOPKG_PACKAGES_HOST_URL="https://global.infrastructureresilience.org/packages"
AUTOPKG_PACKAGES_HOST_URL="https://irv-autopkg.s3.eu-west-2.amazonaws.com"
GDAL_CACHEMAX=2048
AUTOPKG_INCLUDE_TEST_PROCESSORS="True"

# AWSS3 Backend
# Dev
AUTOPKG_TEST_S3_ACCESS_KEY=
AUTOPKG_TEST_S3_SECRET_KEY=
AUTOPKG_TEST_S3_BUCKET=irv-autopkg-dev

# Prod
AUTOPKG_S3_ACCESS_KEY=
AUTOPKG_S3_SECRET_KEY=
AUTOPKG_S3_BUCKET=irv-autopkg
AUTOPKG_S3_REGION=eu-west-2

# AWS OSM / Damages DB
AUTOPKG_OSM_PGHOST=
AUTOPKG_OSM_PORT=
AUTOPKG_OSM_PGDATABASE=
AUTOPKG_OSM_PGUSER=
AUTOPKG_OSM_PGPASSWORD=
```

To pull and run images:

```
docker compose -f docker-compose.yaml up -d
```

`docker compose up --build` will build the necessary images and bring up
containers with environment data provided from the .env file.

`docker ps` to see the containers' state.
