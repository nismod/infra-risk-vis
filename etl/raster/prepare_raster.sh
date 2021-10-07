#!/bin/env bash

set -o errexit

# TODO validate arguments

INPUT_FILE=$1
OUTPUT_FILE=$2

NODATA=$(gdalinfo "$INPUT_FILE" -json | jq .bands[0].noDataValue)

# handle case of NODATA == nan - the JSON output of gdalinfo will change nan to "NaN"
# so we need to reverse that in order for gdal_calc.py to not error
if [ "$NODATA" == '"NaN"' ]; then
  NODATA=nan
fi

# TODO more robust handling of tmp file
# TODO add CLI toggle for removing zeros (for some datasets you might want to leave the 0 values intact)

# replace zeros with NoData value
gdal_calc.py -A "$INPUT_FILE" --outfile=tmp.tif --overwrite --calc="numpy.where(A==0,$NODATA,A)" --hideNoData --NoDataValue=$NODATA

# TODO tune COG conversion parameters

# translate to Cloud-Optimised GeoTIFF
gdal_translate tmp.tif "$OUTPUT_FILE" -of COG

rm tmp.tif