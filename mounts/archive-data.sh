#!/bin/bash

set -ex

DATESTR=$(date -I)
tar cJf "${DATESTR}_tileserver-vector-data.tar.xz" vector/
tar cJf "${DATESTR}_tileserver-raster-data.tar.xz" raster/
