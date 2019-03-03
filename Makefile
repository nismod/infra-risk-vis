.PHONY: all air boundaries rail road water clean

in = incoming_data
temp = intermediate_data
out = data

all: air boundaries rail road water

air: $(out)/air.mbtiles
boundaries: $(out)/boundaries.mbtiles
rail: $(out)/rail.mbtiles
road: $(out)/road.mbtiles
water: $(out)/water.mbtiles

$(out)/air.mbtiles: $(in)/network/air_edges.shp $(in)/network/air_edges.shp
	ogr2ogr -f GeoJSON $(temp)/air_edges.json -t_srs EPSG:4326 $(in)/network/air_edges.shp
	
	rm -f $(out)/air.mbtiles
	tippecanoe \
    --no-feature-limit \
    --no-line-simplification \
    --no-tile-size-limit \
    -o $(out)/air.mbtiles \
    $(temp)/air_edges.json

$(out)/boundaries.mbtiles: $(in)/boundaries/admin_0_boundaries.shp $(in)/boundaries/admin_0_boundaries_labels.shp $(in)/boundaries/admin_1_boundaries.shp $(in)/boundaries/admin_1_boundaries_labels.shp $(in)/boundaries/physical_lakes.shp 
	ogr2ogr -f GeoJSON $(temp)/admin_0_boundaries.json -t_srs EPSG:4326 $(in)/boundaries/admin_0_boundaries.shp
	ogr2ogr -f GeoJSON $(temp)/admin_0_boundaries_labels.json -t_srs EPSG:4326 $(in)/boundaries/admin_0_boundaries_labels.shp
	ogr2ogr -f GeoJSON $(temp)/admin_1_boundaries.json -t_srs EPSG:4326 $(in)/boundaries/admin_1_boundaries.shp
	ogr2ogr -f GeoJSON $(temp)/admin_1_boundaries_labels.json -t_srs EPSG:4326 $(in)/boundaries/admin_1_boundaries_labels.shp
	ogr2ogr -f GeoJSON $(temp)/physical_lakes.json -t_srs EPSG:4326 $(in)/boundaries/physical_lakes.shp
	
	rm -f $(out)/boundaries.mbtiles
	tippecanoe \
	-zg \
	-r1 \
    --no-feature-limit \
    --no-line-simplification \
    --no-tile-size-limit \
	--extend-zooms-if-still-dropping \
    -o $(out)/boundaries.mbtiles \
    $(temp)/admin_0_boundaries.json \
	$(temp)/admin_0_boundaries_labels.json \
    $(temp)/admin_1_boundaries.json \
	$(temp)/admin_1_boundaries_labels.json \
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
    -o $(out)/road.mbtiles \
    $(temp)/road_edges_national.json \
    $(temp)/road_edges_provincial.json \
	$(temp)/road_edges_rural.json

$(out)/water.mbtiles: $(in)/network/water_edges.shp $(in)/network/water_nodes.shp
	ogr2ogr -f GeoJSON $(temp)/water_edges.json -t_srs EPSG:4326 $(in)/network/water_edges.shp
	ogr2ogr -f GeoJSON $(temp)/water_nodes.json -t_srs EPSG:4326 $(in)/network/water_nodes.shp

	rm -f $(out)/water.mbtiles
	tippecanoe \
    -o $(out)/water.mbtiles \
    $(temp)/water_edges.json \
    $(temp)/water_nodes.json \

clean:
	rm -f $(out)/*.mbtiles