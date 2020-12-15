#!/usr/bin/env bash
set -e
set -x

#
# Deploy app
#

# built files for frontend
rsync -rvz build/ raghavpant@tool.oi-analytics.com:/var/www/html

# data and config for tileserver
rsync -rvz data/ raghavpant@tool.oi-analytics.com:/var/www/tileserver/data
rsync -rvz styles/ raghavpant@tool.oi-analytics.com:/var/www/tileserver/styles
rsync -rvz fonts/ raghavpant@tool.oi-analytics.com:/var/www/tileserver/fonts
rsync -rvz config.json raghavpant@tool.oi-analytics.com:/var/www/tileserver
