.PHONY: all clean

all: ./data/pot_edges.mbtiles ./data/abs_nodes.mbtiles ./data/rail_edges.mbtiles ./data/rail_nodes.mbtiles ./data/road_edges.mbtiles ./data/bridges.mbtiles ./data/elec_edges.mbtiles ./data/elec_nodes.mbtiles

./data/rail_edges.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--output=./data/rail_edges.mbtiles \
		--layer=rail_edges \
		./intermediate_data/rail_edges.json

./data/rail_nodes.mbtiles:
	tippecanoe \
		-zg \
		--output=./data/rail_nodes.mbtiles \
		--layer=rail_nodes \
		./intermediate_data/rail_nodes.json

./data/road_edges.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./data/road_edges.mbtiles \
		--layer=road_edges \
		./intermediate_data/road_edges.json

./data/bridges.mbtiles:
	tippecanoe \
		-zg \
		--output=./data/bridges.mbtiles \
		--layer=bridges \
		./intermediate_data/bridges.json

./data/elec_edges.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./data/elec_edges.mbtiles \
		--layer=elec_edges \
		./intermediate_data/elec_edges.json

./data/elec_nodes.mbtiles:
	tippecanoe \
		-zg \
		--output=./data/elec_nodes.mbtiles \
		--layer=elec_nodes \
		./intermediate_data/elec_nodes.json

./data/pot_edges.mbtiles:
	tippecanoe \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./data/pot_edges.mbtiles \
		--layer=pot_edges \
		./intermediate_data/pot_edges.json

./data/abs_nodes.mbtiles:
	tippecanoe \
		-zg \
		--output=./data/abs_nodes.mbtiles \
		--layer=abs_nodes \
		./intermediate_data/abs_nodes.json

clean:
	rm -f ./data/*.mbtiles
