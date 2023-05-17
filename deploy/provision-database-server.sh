#!/usr/bin/env bash

#
# Provision virtual machine as database server
# - assuming OS is Ubuntu 20.04 LTS
#

# Configure postgres and postgis versions
export IMAGE_VERSION=$(lsb_release -cs)
export POSTGRES_MAJOR_VERSION=14
export POSTGIS_MAJOR_VERSION=3

# Install helper packages
sudo apt-get -qq update --fix-missing && sudo apt-get -qq --yes upgrade

export DEBIAN_FRONTEND=noninteractive \
    && sudo apt-get update \
    && sudo apt-get -y install \
        gnupg2 wget ca-certificates pwgen software-properties-common \
        apt-transport-https curl gettext

# Install postgres
sudo sh -c "echo \"deb http://apt.postgresql.org/pub/repos/apt/ ${IMAGE_VERSION}-pgdg main\" > /etc/apt/sources.list.d/pgdg.list"
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc -O- | sudo apt-key add -


export DEBIAN_FRONTEND=noninteractive \
    && sudo apt-get update \
    && sudo apt-get -y --no-install-recommends install \
      postgresql-client-${POSTGRES_MAJOR_VERSION} \
        postgresql-common \
        postgresql-${POSTGRES_MAJOR_VERSION} \
        postgresql-${POSTGRES_MAJOR_VERSION}-postgis-${POSTGIS_MAJOR_VERSION} \
        postgresql-${POSTGRES_MAJOR_VERSION}-postgis-${POSTGIS_MAJOR_VERSION}-scripts

# Create database
sudo -u postgres createdb jamaicadev
# Create admin user for create/backup/restore
sudo -u postgres createuser jamaicadev --pwprompt
sudo -u postgres psql -c 'alter user jamaicadev with superuser;'
# Create read-only user for app connection
sudo -u postgres createuser jsratapp --pwprompt
sudo -u postgres psql -c 'GRANT pg_read_all_data TO jsratapp;''
