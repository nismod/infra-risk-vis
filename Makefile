.PHONY: all clean

all: ./tileserver/data/pot_edges.mbtiles ./tileserver/data/abs_nodes.mbtiles ./tileserver/data/rail_edges.mbtiles ./tileserver/data/rail_nodes.mbtiles ./tileserver/data/road_edges.mbtiles ./tileserver/data/bridges.mbtiles ./tileserver/data/elec_edges.mbtiles ./tileserver/data/elec_nodes.mbtiles

./tileserver/data/rail_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--output=./tileserver/data/rail_edges.mbtiles \
		--layer=rail_edges \
		./intermediate_data/rail_edges.json

./tileserver/data/rail_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=./tileserver/data/rail_nodes.mbtiles \
		--layer=rail_nodes \
		./intermediate_data/rail_nodes.json

./tileserver/data/road_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./tileserver/data/road_edges.mbtiles \
		--layer=road_edges \
		./intermediate_data/road_edges.json

./tileserver/data/bridges.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=./tileserver/data/bridges.mbtiles \
		--layer=bridges \
		./intermediate_data/bridges.json

./tileserver/data/elec_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./tileserver/data/elec_edges.mbtiles \
		--layer=elec_edges \
		./intermediate_data/elec_edges.json

./tileserver/data/elec_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=./tileserver/data/elec_nodes.mbtiles \
		--layer=elec_nodes \
		./intermediate_data/elec_nodes.json

./tileserver/data/pot_edges.mbtiles:
	tippecanoe \
		--generate-ids \
		--minimum-zoom=3 \
		--maximum-zoom=15 \
		--drop-smallest-as-needed \
		--output=./tileserver/data/pot_edges.mbtiles \
		--layer=pot_edges \
		./intermediate_data/pot_edges.json

./tileserver/data/abs_nodes.mbtiles:
	tippecanoe \
		-zg \
		--generate-ids \
		--output=./tileserver/data/abs_nodes.mbtiles \
		--layer=abs_nodes \
		./intermediate_data/abs_nodes.json

clean:
	rm -f ./tileserver/data/*.mbtiles
