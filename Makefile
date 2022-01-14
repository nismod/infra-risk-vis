# load variables from .env file
ifneq (,$(wildcard ./.env))
    include .env
    export
endif

# get path to this Makefile regardless of CWD - see https://stackoverflow.com/a/18137056/1478817
current_dir := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))


.PHONY: all clean vector networks raster raster-fluvial raster-surface raster-coastal raster-cyclone clean-vector clean-rasters boundaries boundaries-parish boundaries-community

all: vector raster


vector: networks boundaries

networks:
	"$(current_dir)/etl/vector/process_all.py" --raw "${RAW_DATA_DIR}" --out "$(current_dir)/tileserver/vector/data"

clean-vector:
	rm ./tileserver/data/*.mbtiles


boundaries-parish:
	"$(current_dir)/etl/vector/process_boundaries.sh" "${RAW_DATA_DIR}/boundaries/admin_boundaries.gpkg" "$(current_dir)/tileserver/vector/data/boundaries_parish.mbtiles" parish admin1

boundaries-community:
	"$(current_dir)/etl/vector/process_boundaries.sh" "${RAW_DATA_DIR}/boundaries/admin_boundaries.gpkg" "$(current_dir)/tileserver/vector/data/boundaries_community.mbtiles" community admin2

boundaries: boundaries-parish boundaries-community


RASTER_BASE_COMMAND = "$(current_dir)/etl/raster/process_all.py" --raw "${RAW_DATA_DIR}" --out "$(current_dir)/tileserver/raster/input"

raster:
	${RASTER_BASE_COMMAND}

raster-fluvial:
	${RASTER_BASE_COMMAND} --type fluvial

raster-surface:
	${RASTER_BASE_COMMAND} --type surface

raster-coastal:
	${RASTER_BASE_COMMAND} --type coastal

raster-cyclone:
	${RASTER_BASE_COMMAND} --type cyclone


clean-raster:
	rm ./tileserver-raster/data/*.tif


clean: clean-vector clean-raster
