.PHONY: all clean vector raster clean-vector clean-rasters

all: vector raster

vector_layers = boundaries rail_edges rail_nodes road_edges bridges elec_edges elec_nodes pot_edges abs_nodes
raster_layers = fluvial_rp20_raw fluvial_rp50_raw fluvial_rp100_raw fluvial_rp200_raw fluvial_rp500_raw fluvial_rp1500_raw coastal_rp1_raw coastal_rp2_raw coastal_rp5_raw coastal_rp10_raw coastal_rp50_raw coastal_rp100_raw

vector: $(patsubst %,./tileserver/data/%.mbtiles,$(vector_layers))
clean-vector:
	rm ./tileserver/data/*.mbtiles

raster: $(patsubst %,./tileserver-raster/data/%.tif,$(raster_layers))
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

PROCESS_RASTER = ./scripts/raster/prepare_raster.sh

./tileserver-raster/data/fluvial_rp20_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q20_RD_02.tif $@

./tileserver-raster/data/fluvial_rp50_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q50_RD_02.tif $@

./tileserver-raster/data/fluvial_rp100_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q100_RD_02.tif $@

./tileserver-raster/data/fluvial_rp200_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q200_RD_02.tif $@

./tileserver-raster/data/fluvial_rp500_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q500_RD_02.tif $@

./tileserver-raster/data/fluvial_rp1500_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Fluvial/JM_FLRF_UD_Q1500_RD_02.tif $@


./tileserver-raster/data/coastal_rp1_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_1.tif $@

./tileserver-raster/data/coastal_rp2_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_2.tif $@

./tileserver-raster/data/coastal_rp5_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_5.tif $@

./tileserver-raster/data/coastal_rp10_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_10.tif $@

./tileserver-raster/data/coastal_rp50_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_50.tif $@

./tileserver-raster/data/coastal_rp100_raw.tif:
	$(PROCESS_RASTER) ./intermediate_data/Coastal/JamaicaJAM001RCP452010_epsg_32618_RP_100.tif $@

clean: clean-vector clean-raster
