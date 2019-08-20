#
# Build and deploy app
#

npm run build
rsync -ravz build/ raghav@argentina.oi-analytics.com:/var/www/html
rsync -ravz data/ raghav@argentina.oi-analytics.com:/var/www/tileserver/data
rsync -ravz styles/ raghav@argentina.oi-analytics.com:/var/www/tileserver/styles
