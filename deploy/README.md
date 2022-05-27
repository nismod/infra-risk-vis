# Deploy

The site runs on a single Linux virtual machine using docker compose.

The stack consists of five services:
- Web server (nginx)
- Raster tileserver (terracotta)
- Vector tileserver (tileserver-gl)
- Backend / API (bespoke Python app)
- Database (PostgreSQL with PostGIS)

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

Update the .env file with the relevant database connection parameters:
```
PGDATABASE=jamaica
PGUSER=docker
PGPASSWORD=docker
PGHOST=localhost
```

`docker compose up --build` will build the necessary images and bring up
containers with environment data provided from the .env file.

`docker ps` to see the containers' state.
