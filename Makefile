.PHONY: all clean vector raster

all: vector raster

vector_layers = boundaries rail_edges rail_nodes road_edges bridges elec_edges elec_nodes pot_edges abs_nodes
raster_layers = fluvial_rp20_raw fluvial_rp50_raw fluvial_rp100_raw fluvial_rp200_raw fluvial_rp500_raw fluvial_rp1500_raw

vector: $(patsubst %,./tileserver/data/%.mbtiles,$(vector_layers))
raster: $(patsubst %,./tileserver-raster/data/%.tif,$(raster_layers))

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

MAKE_COG = ./scripts/raster/make_cog.sh

./tileserver-raster/data/fluvial_rp20_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q20_RD_02.tif $@

./tileserver-raster/data/fluvial_rp50_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q50_RD_02.tif $@

./tileserver-raster/data/fluvial_rp100_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q100_RD_02.tif $@

./tileserver-raster/data/fluvial_rp200_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q200_RD_02.tif $@

./tileserver-raster/data/fluvial_rp500_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q500_RD_02.tif $@

./tileserver-raster/data/fluvial_rp1500_raw.tif:
	$(MAKE_COG) ./intermediate_data/Fluvial/JM_FLRF_UD_Q1500_RD_02.tif $@

clean:
	rm -f ./tileserver/data/*.mbtiles && rm -f ./tileserver-raster/data/*.tif
