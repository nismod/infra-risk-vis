.PHONY: all air boundaries bridge flood rail road water clean

all: flood road boundaries

air: ./data/air_edges.mbtiles ./data/air_nodes.mbtiles
boundaries: ./data/boundaries.mbtiles
rail: ./data/rail_edges.mbtiles ./data/rail_nodes.mbtiles
road: ./data/road_edges.mbtiles
bridge: ./data/bridges.mbtiles
water: ./data/port_edges.mbtiles ./data/port_nodes.mbtiles
flood: ./data/flood.mbtiles

./data/air_edges.mbtiles: ./intermediate_data/air_edges.json ./intermediate_data/air_edges.json
	rm -f ./data/air_edges.mbtiles
	tippecanoe \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		-o ./data/air_edges.mbtiles \
		./intermediate_data/air_edges.json

./data/air_nodes.mbtiles: ./intermediate_data/air_nodes.json ./intermediate_data/air_nodes.json
	rm -f ./data/air_nodes.mbtiles
	tippecanoe \
		-zg \
		-r0 \
		--no-feature-limit \
		--no-tile-size-limit \
		-o ./data/air_nodes.mbtiles \
		./intermediate_data/air_nodes.json


./data/boundaries.mbtiles: ./intermediate_data/boundaries/CNTR_LB_2020_4326.json ./intermediate_data/boundaries/WB_Admin0_disputed_areas.geojson ./intermediate_data/boundaries/WB_countries_Admin0.geojson
	rm -f ./data/boundaries.mbtiles
	ogr2ogr -f GeoJSON \
		./intermediate_data/boundaries/WB_Admin0_disputed_areas.json \
		-t_srs EPSG:4326 \
		./intermediate_data/boundaries/WB_Admin0_disputed_areas.geojson
	ogr2ogr -f GeoJSON \
		./intermediate_data/boundaries/WB_countries_Admin0.json \
		-t_srs EPSG:4326 \
		./intermediate_data/boundaries/WB_countries_Admin0.geojson
	tippecanoe \
		-zg \
		-r1 -pk -pf \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o ./data/boundaries.mbtiles \
		./intermediate_data/boundaries/WB_Admin0_disputed_areas.json \
		./intermediate_data/boundaries/CNTR_LB_2020_4326.json \
		./intermediate_data/boundaries/WB_countries_Admin0.json

./data/bridges.mbtiles: ./intermediate_data/bridges.json ./intermediate_data/bridges.json
	rm -f ./data/bridges.mbtiles
	tippecanoe \
		-zg \
		-r0 \
		--no-feature-limit \
		--no-tile-size-limit \
		-o ./data/bridges.mbtiles \
		./intermediate_data/bridges.json

./data/rail_nodes.mbtiles: ./intermediate_data/rail_nodes.json ./intermediate_data/rail_nodes.json
	rm -f ./data/rail_nodes.mbtiles
	tippecanoe \
		-o ./data/rail_nodes.mbtiles \
		./intermediate_data/rail_nodes.json

./data/rail_edges.mbtiles: ./intermediate_data/rail_edges.json ./intermediate_data/rail_edges.json
	rm -f ./data/rail_edges.mbtiles
	tippecanoe \
		-o ./data/rail_edges.mbtiles \
		./intermediate_data/rail_edges.json

./data/road_edges.mbtiles: ./intermediate_data/road_edges.json ./intermediate_data/road_edges.json
	rm -f ./data/road_edges.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o ./data/road_edges.mbtiles \
		./intermediate_data/road_edges.json \

./data/port_nodes.mbtiles: ./intermediate_data/port_nodes.json ./intermediate_data/port_nodes.json
	rm -f ./data/port_nodes.mbtiles
	tippecanoe \
		-zg \
		-r0 \
		--no-feature-limit \
		--no-tile-size-limit \
		-o ./data/port_nodes.mbtiles \
		./intermediate_data/port_nodes.json

./data/port_edges.mbtiles: ./intermediate_data/port_edges.json ./intermediate_data/port_edges.json
	rm -f ./data/port_edges.mbtiles
	tippecanoe \
		-o ./data/port_edges.mbtiles \
		./intermediate_data/port_edges.json

./data/flood.mbtiles: ./intermediate_data/flood_data/baseline_fluvial_1in1000_1m2m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_2m3m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_3m4m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_4m999m.json ./intermediate_data/flood_data/rcp_4p5_fluvial_1in1000_1m2m.json ./intermediate_data/flood_data/rcp_4p5_fluvial_1in1000_2m3m.json ./intermediate_data/flood_data/rcp_4p5_fluvial_1in1000_3m4m.json ./intermediate_data/flood_data/rcp_4p5_fluvial_1in1000_4m999m.json ./intermediate_data/flood_data/rcp_8p5_fluvial_1in1000_1m2m.json ./intermediate_data/flood_data/rcp_8p5_fluvial_1in1000_2m3m.json ./intermediate_data/flood_data/rcp_8p5_fluvial_1in1000_3m4m.json ./intermediate_data/flood_data/rcp_8p5_fluvial_1in1000_4m999m.json 
	rm -f ./data/flood.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o ./data/flood.mbtiles \
		./intermediate_data/flood_data/*

clean:
	rm -f ./data/*.mbtiles
