# load variables from .env file
ifneq (,$(wildcard ./.env))
    include .env
    export
endif

# get path to this Makefile regardless of CWD - see https://stackoverflow.com/a/18137056/1478817
current_dir := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))


.PHONY: all clean vector networks raster raster-fluvial raster-surface raster-coastal raster-cyclone clean-vector clean-rasters boundaries boundaries-parish boundaries-community boundaries-subdivision

all: vector raster


# ======
# Networks vector data

vector: networks boundaries

networks:
	"$(current_dir)/etl/vector/process_all.py" --raw "${RAW_DATA_DIR}" --out "$(current_dir)/tileserver/vector/data"

clean-vector:
	rm ./tileserver/data/*.mbtiles


# ======
# Administrative boundaries

BOUNDARY_COMMAND = "$(current_dir)/etl/vector/process_boundaries.sh"
BOUNDARY_FILE = "${RAW_DATA_DIR}/boundaries/admin_boundaries.gpkg"
BOUNDARY_OUT_BASE_PATH = "$(current_dir)/tileserver/vector/data"

boundaries-parish:
	 "$(BOUNDARY_COMMAND)" "$(BOUNDARY_FILE)" "$(BOUNDARY_OUT_BASE_PATH)/boundaries_parish.mbtiles" parish admin1

boundaries-community:
	"$(BOUNDARY_COMMAND)" "$(BOUNDARY_FILE)" "$(BOUNDARY_OUT_BASE_PATH)/boundaries_community.mbtiles" community admin2

boundaries-subdivision:
	"$(BOUNDARY_COMMAND)" "$(BOUNDARY_FILE)" "$(BOUNDARY_OUT_BASE_PATH)/boundaries_subdivision.mbtiles" community admin3

boundaries: boundaries-parish boundaries-community boundaries-subdivision


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


clean: clean-vector clean-raster
