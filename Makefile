.PHONY: all air boundaries bridge flood rail road water clean

all: air boundaries bridge flood rail road water

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
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		-o ./data/air_nodes.mbtiles \
		./intermediate_data/air_nodes.json

./incoming_data/boundaries/admin_0_labels.shp:  ./incoming_data/boundaries/admin_0_labels.csv
	ogr2ogr -f "ESRI Shapefile" $@ $< -a_srs EPSG:4326 \
		-oo X_POSSIBLE_NAMES=lon \
		-oo Y_POSSIBLE_NAMES=lat \
		-oo KEEP_GEOM_COLUMNS=NO

./incoming_data/boundaries/admin_1_labels.shp:  ./incoming_data/boundaries/admin_1_labels.csv
	ogr2ogr -f "ESRI Shapefile" $@ $< -a_srs EPSG:4326 \
		-oo X_POSSIBLE_NAMES=lon \
		-oo Y_POSSIBLE_NAMES=lat \
		-oo KEEP_GEOM_COLUMNS=NO

./data/boundaries.mbtiles: ./incoming_data/boundaries/admin_0_boundaries.shp ./incoming_data/boundaries/admin_0_labels.shp ./incoming_data/boundaries/admin_1_boundaries.shp ./incoming_data/boundaries/admin_1_labels.shp ./incoming_data/boundaries/physical_lakes.shp
	ogr2ogr -f GeoJSON \
		./intermediate_data/admin_0_boundaries.json \
		./incoming_data/boundaries/admin_0_boundaries.shp
	ogr2ogr -f GeoJSON \
		./intermediate_data/admin_0_labels.json \
		-t_srs EPSG:4326 \
		./incoming_data/boundaries/admin_0_labels.shp
	ogr2ogr -f GeoJSON \
		./intermediate_data/admin_1_boundaries.json \
		-t_srs EPSG:4326 \
		./incoming_data/boundaries/admin_1_boundaries.shp
	ogr2ogr -f GeoJSON \
		./intermediate_data/admin_1_labels.json \
		-t_srs EPSG:4326 \
		./incoming_data/boundaries/admin_1_labels.shp
	ogr2ogr -f GeoJSON \
		./intermediate_data/physical_lakes.json \
		-t_srs EPSG:4326 \
		./incoming_data/boundaries/physical_lakes.shp
	rm -f ./data/boundaries.mbtiles
	tippecanoe \
		-zg \
		-r1 -pk -pf \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o ./data/boundaries.mbtiles \
		./intermediate_data/admin_0_boundaries.json \
		./intermediate_data/admin_0_labels.json \
		./intermediate_data/admin_1_boundaries.json \
		./intermediate_data/admin_1_labels.json \
		./intermediate_data/physical_lakes.json

./data/bridges.mbtiles: ./intermediate_data/bridges.json ./intermediate_data/bridges.json
	rm -f ./data/bridges.mbtiles
	tippecanoe \
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
		-o ./data/port_nodes.mbtiles \
		./intermediate_data/port_nodes.json

./data/port_edges.mbtiles: ./intermediate_data/port_edges.json ./intermediate_data/port_edges.json
	rm -f ./data/port_edges.mbtiles
	tippecanoe \
		-o ./data/port_edges.mbtiles \
		./intermediate_data/port_edges.json

./intermediate_data/flood_data/baseline_fluvial_1in1000_1m-2m.json: ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_fluvial_1in1000_1m-2m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_fluvial_1in1000_2m-3m.json: ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_fluvial_1in1000_2m-3m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_fluvial_1in1000_3m-4m.json: ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_fluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_fluvial_1in1000_4m-999m.json: ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_fluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


./intermediate_data/flood_data/baseline_pluvial_1in1000_1m-2m.json: ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_pluvial_1in1000_1m-2m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_pluvial_1in1000_2m-3m.json: ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_pluvial_1in1000_2m-3m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_pluvial_1in1000_3m-4m.json: ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_pluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/baseline_pluvial_1in1000_4m-999m.json: ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/baseline_pluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Baseline/pluvial/P_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


./intermediate_data/flood_data/high_fluvial_1in1000_1m-2m.json: ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_fluvial_1in1000_1m-2m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_fluvial_1in1000_2m-3m.json: ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_fluvial_1in1000_2m-3m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_fluvial_1in1000_3m-4m.json: ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_fluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_fluvial_1in1000_4m-999m.json: ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_fluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


./intermediate_data/flood_data/high_pluvial_1in1000_1m-2m.json: ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_pluvial_1in1000_1m-2m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_pluvial_1in1000_2m-3m.json: ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_pluvial_1in1000_2m-3m.json -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_pluvial_1in1000_3m-4m.json: ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_pluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

./intermediate_data/flood_data/high_pluvial_1in1000_4m-999m.json: ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON ./intermediate_data/flood_data/high_pluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 ./incoming_data/flood_data/FATHOM/Future_High/pluvial/P_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


./data/flood.mbtiles: ./intermediate_data/flood_data/baseline_fluvial_1in1000_1m-2m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_2m-3m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_3m-4m.json ./intermediate_data/flood_data/baseline_fluvial_1in1000_4m-999m.json ./intermediate_data/flood_data/baseline_pluvial_1in1000_1m-2m.json ./intermediate_data/flood_data/baseline_pluvial_1in1000_2m-3m.json ./intermediate_data/flood_data/baseline_pluvial_1in1000_3m-4m.json ./intermediate_data/flood_data/baseline_pluvial_1in1000_4m-999m.json ./intermediate_data/flood_data/high_fluvial_1in1000_1m-2m.json ./intermediate_data/flood_data/high_fluvial_1in1000_2m-3m.json ./intermediate_data/flood_data/high_fluvial_1in1000_3m-4m.json ./intermediate_data/flood_data/high_fluvial_1in1000_4m-999m.json ./intermediate_data/flood_data/high_pluvial_1in1000_1m-2m.json ./intermediate_data/flood_data/high_pluvial_1in1000_2m-3m.json ./intermediate_data/flood_data/high_pluvial_1in1000_3m-4m.json ./intermediate_data/flood_data/high_pluvial_1in1000_4m-999m.json
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
