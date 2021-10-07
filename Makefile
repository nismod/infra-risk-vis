# load variables from .env file
ifneq (,$(wildcard ./.env))
    include .env
    export
endif

# get path to this Makefile regardless of CWD - see https://stackoverflow.com/a/18137056/1478817
current_dir := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
# current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))


.PHONY: all clean vector raster raster-fluvial raster-surface raster-coastal raster-cyclone clean-vector clean-rasters

all: vector raster

vector_layers = boundaries rail_edges rail_nodes road_edges bridges elec_edges elec_nodes pot_edges abs_nodes
raster_layers = fluvial_rp20_raw fluvial_rp50_raw fluvial_rp100_raw fluvial_rp200_raw fluvial_rp500_raw fluvial_rp1500_raw coastal_rp1_raw coastal_rp2_raw coastal_rp5_raw coastal_rp10_raw coastal_rp50_raw coastal_rp100_raw

vector: $(patsubst %,./tileserver/data/%.mbtiles,$(vector_layers))
clean-vector:
	rm ./tileserver/data/*.mbtiles

RASTER_BASE_COMMAND = "$(current_dir)/etl/raster/process_all.py" --raw "${RAW_DATA_DIR}" --out "$(current_dir)/tileserver/raster/data"

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

./tileserver/data/boundaries.mbtiles:
	cp ./incoming_data/boundaries.mbtiles $@

./tileserver/data/rail_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--output=$@ \
		--layer=rail_edges \
		./intermediate_data/rail_edges.json

./tileserver/data/rail_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=$@ \
		--layer=rail_nodes \
		./intermediate_data/rail_nodes.json

./tileserver/data/road_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=$@ \
		--layer=road_edges \
		./intermediate_data/road_edges.json

./tileserver/data/bridges.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=$@ \
		--layer=bridges \
		./intermediate_data/bridges.json

./tileserver/data/elec_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=$@ \
		--layer=elec_edges \
		./intermediate_data/elec_edges.json

./tileserver/data/elec_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=$@ \
		--layer=elec_nodes \
		./intermediate_data/elec_nodes.json

./tileserver/data/pot_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=$@ \
		--layer=pot_edges \
		./intermediate_data/pot_edges.json

./tileserver/data/abs_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=$@ \
		--layer=abs_nodes \
		./intermediate_data/abs_nodes.json


clean: clean-vector clean-raster
