.PHONY: all air boundaries rail road water flood clean

in = incoming_data
temp = intermediate_data
out = data

all: air boundaries rail road water flood

air: $(out)/air.mbtiles
boundaries: $(out)/boundaries.mbtiles
rail: $(out)/rail.mbtiles
road: $(out)/road.mbtiles
water: $(out)/water.mbtiles
flood: $(out)/flood.mbtiles

$(out)/air.mbtiles: $(in)/network/air_edges.shp $(in)/network/air_edges.shp
	ogr2ogr -f GeoJSON $(temp)/air_edges.json -t_srs EPSG:4326 $(in)/network/air_edges.shp
	rm -f $(out)/air.mbtiles
	tippecanoe \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		-o $(out)/air.mbtiles \
		$(temp)/air_edges.json

 $(in)/boundaries/admin_0_labels.shp:  $(in)/boundaries/admin_0_labels.csv
	ogr2ogr -f "ESRI Shapefile" $@ $< -a_srs EPSG:4326 \
		-oo X_POSSIBLE_NAMES=lon \
		-oo Y_POSSIBLE_NAMES=lat \
		-oo KEEP_GEOM_COLUMNS=NO

 $(in)/boundaries/admin_1_labels.shp:  $(in)/boundaries/admin_1_labels.csv
	ogr2ogr -f "ESRI Shapefile" $@ $< -a_srs EPSG:4326 \
		-oo X_POSSIBLE_NAMES=lon \
		-oo Y_POSSIBLE_NAMES=lat \
		-oo KEEP_GEOM_COLUMNS=NO

$(out)/boundaries.mbtiles: $(in)/boundaries/admin_0_boundaries.shp $(in)/boundaries/admin_0_labels.shp $(in)/boundaries/admin_1_boundaries.shp $(in)/boundaries/admin_1_labels.shp $(in)/boundaries/physical_lakes.shp
	ogr2ogr -f GeoJSON \
		$(temp)/admin_0_boundaries.json \
		$(in)/boundaries/admin_0_boundaries.shp
	ogr2ogr -f GeoJSON \
		$(temp)/admin_0_labels.json \
		-t_srs EPSG:4326 \
		$(in)/boundaries/admin_0_labels.shp
	ogr2ogr -f GeoJSON \
		$(temp)/admin_1_boundaries.json \
		-t_srs EPSG:4326 \
		$(in)/boundaries/admin_1_boundaries.shp
	ogr2ogr -f GeoJSON \
		$(temp)/admin_1_labels.json \
		-t_srs EPSG:4326 \
		$(in)/boundaries/admin_1_labels.shp
	ogr2ogr -f GeoJSON \
		$(temp)/physical_lakes.json \
		-t_srs EPSG:4326 \
		$(in)/boundaries/physical_lakes.shp
	rm -f $(out)/boundaries.mbtiles
	tippecanoe \
		-zg \
		-r1 -pk -pf \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o $(out)/boundaries.mbtiles \
		$(temp)/admin_0_boundaries.json \
		$(temp)/admin_0_labels.json \
		$(temp)/admin_1_boundaries.json \
		$(temp)/admin_1_labels.json \
		$(temp)/physical_lakes.json

$(out)/rail.mbtiles: $(in)/network/rail_edges.shp $(in)/network/rail_edges.shp
	ogr2ogr -f GeoJSON $(temp)/rail_edges.json -t_srs EPSG:4326 $(in)/network/rail_edges.shp

	rm -f $(out)/rail.mbtiles
	tippecanoe \
		-o $(out)/rail.mbtiles \
		$(temp)/rail_edges.json

$(out)/road.mbtiles: $(in)/network/road_edges_national.shp $(in)/network/road_edges_provincial.shp $(in)/network/road_edges_rural.shp
	ogr2ogr -f GeoJSON $(temp)/road_edges_national.json -t_srs EPSG:4326 $(in)/network/road_edges_national.shp
	ogr2ogr -f GeoJSON $(temp)/road_edges_provincial.json -t_srs EPSG:4326 $(in)/network/road_edges_provincial.shp
	ogr2ogr -f GeoJSON $(temp)/road_edges_rural.json -t_srs EPSG:4326 $(in)/network/road_edges_rural.shp

	rm -f $(out)/road.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o $(out)/road.mbtiles \
		$(temp)/road_edges_national.json \
		$(temp)/road_edges_provincial.json \
		$(temp)/road_edges_rural.json

$(out)/water.mbtiles: $(in)/network/water_edges.shp $(in)/network/water_nodes.shp
	ogr2ogr -f GeoJSON $(temp)/water_edges.json -t_srs EPSG:4326 $(in)/network/water_edges.shp -s_srs EPSG:4326
	ogr2ogr -f GeoJSON $(temp)/water_nodes.json -t_srs EPSG:4326 $(in)/network/water_nodes.shp -s_srs EPSG:4326

	rm -f $(out)/water.mbtiles
	tippecanoe \
		-o $(out)/water.mbtiles \
		$(temp)/water_edges.json \
		$(temp)/water_nodes.json \

FLOOD_DIRS = $(shell find $(in)/flood_data/ -type d)
FLOOD_FILES = $(shell find $(in)/flood_data/ -type f -name '*')
$(out)/flood.mbtiles: $(in)/flood_data/ $(FLOOD_DIRS) $(FLOOD_FILES)
	mkdir -p $(temp)/flood_data/FATHOM/Baseline/fluvial/
	ogr2ogr -f GeoJSON $(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_1m-2m_threshold.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_1m-2m_threshold.shp -s_srs EPSG:4326
	ogr2ogr -f GeoJSON $(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_2m-3m_threshold.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_2m-3m_threshold.shp -s_srs EPSG:4326
	ogr2ogr -f GeoJSON $(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_3m-4m_threshold.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_3m-4m_threshold.shp -s_srs EPSG:4326
	ogr2ogr -f GeoJSON $(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_4m-999m_threshold.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_4m-999m_threshold.shp -s_srs EPSG:4326

	rm -f $(out)/flood.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o $(out)/flood.mbtiles \
		$(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_1m-2m_threshold.json \
		$(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_2m-3m_threshold.json \
		$(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_3m-4m_threshold.json \
		$(temp)/flood_data/FATHOM/Baseline/fluvial/FU_1in500_4m-999m_threshold.json

clean:
	rm -f $(out)/*.mbtiles
