#!/usr/bin/env bash
set -e
set -x

#
# Deploy app
#

# built files for frontend
rsync -ravz build/ raghav@argentina.oi-analytics.com:/var/www/html

# data and config for tileserver
rsync -ravz data/ raghav@argentina.oi-analytics.com:/var/www/tileserver/data
rsync -ravz styles/ raghav@argentina.oi-analytics.com:/var/www/tileserver/styles
rsync -ravz fonts/ raghav@argentina.oi-analytics.com:/var/www/tileserver/fonts
rsync -ravz config.json raghav@argentina.oi-analytics.com:/var/www/tileserver
