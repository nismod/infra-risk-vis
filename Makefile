# load variables from .env file
ifneq (,$(wildcard ./.env))
    include .env
    export
endif

# get path to this Makefile regardless of CWD - see https://stackoverflow.com/a/18137056/1478817
current_dir := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))


.PHONY: all clean networks raster regions raster-fluvial raster-surface raster-coastal raster-cyclone clean-networks clean-rasters regions-parish regions-enumeration

all: networks raster regions


# ======
# Networks vector data

networks:
	"$(current_dir)/etl/vector/process_all.py" --raw "${RAW_DATA_DIR}" --out "$(current_dir)/tileserver/vector/data"

clean-networks:
	rm ./tileserver/data/*.mbtiles


# ======
# Administrative regions

REGION_COMMAND = "$(current_dir)/etl/vector/process_regions.sh"
REGION_OUT_BASE_PATH = "$(current_dir)/tileserver/vector/data"

regions-parish:
	 "$(REGION_COMMAND)" "${RAW_DATA_DIR}/regions/regions_parish.gpkg" "$(REGION_OUT_BASE_PATH)/regions_parish.mbtiles" "$(REGION_OUT_BASE_PATH)/regions_parish_labels.mbtiles" polygons

regions-enumeration:
	"$(REGION_COMMAND)" "${RAW_DATA_DIR}/regions/regions_enumeration.gpkg" "$(REGION_OUT_BASE_PATH)/regions_enumeration.mbtiles" "$(REGION_OUT_BASE_PATH)/regions_enumeration_labels.mbtiles" polygons

regions: regions-parish regions-enumeration


# ======
# Hazard raster data

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


clean: clean-networks clean-raster
