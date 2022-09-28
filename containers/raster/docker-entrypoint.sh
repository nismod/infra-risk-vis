#!/bin/bash

if [ "$1" = 'serve' ]; then

  if [ ! -f "$TC_DRIVER_PATH" ]; then
    echo "Cannot find DB at, aborting: " $TC_DRIVER_PATH
    exit -1
  fi
  if [ $DEPLOYMENT_ENV = "prod" ]; then
    gunicorn --workers 2 --bind 0.0.0.0:$TC_PORT --umask 007 --access-logfile '-' terracotta.server.app:app
  else
    terracotta -c /opt/config.toml serve -d $TC_DRIVER_PATH --allow-all-ips --port $TC_PORT
  fi
  exit 0
fi

exec "$@"
