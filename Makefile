.PHONY: all air boundaries rail road water flood clean

in = incoming_data
temp = intermediate_data
out = data

all: air boundaries rail road water flood bridge

air: $(out)/air_edges.mbtiles $(out)/air_nodes.mbtiles
boundaries: $(out)/boundaries.mbtiles
rail: $(out)/rail_edges.mbtiles $(out)/rail_nodes.mbtiles
road: $(out)/road_edges.mbtiles
bridge: $(out)/bridges.mbtiles
water: $(out)/port_edges.mbtiles $(out)/port_nodes.mbtiles
flood: $(out)/flood.mbtiles

$(out)/air_edges.mbtiles: $(temp)/air_edges.json $(temp)/air_edges.json
	rm -f $(out)/air_edges.mbtiles
	tippecanoe \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		-o $(out)/air_edges.mbtiles \
		$(temp)/air_edges.json

$(out)/air_nodes.mbtiles: $(temp)/air_nodes.json $(temp)/air_nodes.json
	rm -f $(out)/air_nodes.mbtiles
	tippecanoe \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		-o $(out)/air_nodes.mbtiles \
		$(temp)/air_nodes.json

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

$(out)/bridges.mbtiles: $(temp)/bridges.json $(temp)/bridges.json
	rm -f $(out)/bridges.mbtiles
	tippecanoe \
		-o $(out)/bridges.mbtiles \
		$(temp)/bridges.json

$(out)/rail_nodes.mbtiles: $(temp)/rail_nodes.json $(temp)/rail_nodes.json
	rm -f $(out)/rail_nodes.mbtiles
	tippecanoe \
		-o $(out)/rail_nodes.mbtiles \
		$(temp)/rail_nodes.json

$(out)/rail_edges.mbtiles: $(temp)/rail_edges.json $(temp)/rail_edges.json
	rm -f $(out)/rail_edges.mbtiles
	tippecanoe \
		-o $(out)/rail_edges.mbtiles \
		$(temp)/rail_edges.json

$(out)/road_edges.mbtiles: $(temp)/road_edges.json $(temp)/road_edges.json
	rm -f $(out)/road_edges.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o $(out)/road_edges.mbtiles \
		$(temp)/road_edges.json \

$(out)/port_nodes.mbtiles: $(temp)/port_nodes.json $(temp)/port_nodes.json
	rm -f $(out)/port_nodes.mbtiles
	tippecanoe \
		-o $(out)/port_nodes.mbtiles \
		$(temp)/port_nodes.json

$(out)/port_edges.mbtiles: $(temp)/port_edges.json $(temp)/port_edges.json
	rm -f $(out)/port_edges.mbtiles
	tippecanoe \
		-o $(out)/port_edges.mbtiles \
		$(temp)/port_edges.json

$(temp)/flood_data/baseline_fluvial_1in1000_1m-2m.json: $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_fluvial_1in1000_1m-2m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_fluvial_1in1000_2m-3m.json: $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_fluvial_1in1000_2m-3m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_fluvial_1in1000_3m-4m.json: $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_fluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_fluvial_1in1000_4m-999m.json: $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_fluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/fluvial/FU_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


$(temp)/flood_data/baseline_pluvial_1in1000_1m-2m.json: $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_pluvial_1in1000_1m-2m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_pluvial_1in1000_2m-3m.json: $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_pluvial_1in1000_2m-3m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_pluvial_1in1000_3m-4m.json: $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_pluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/baseline_pluvial_1in1000_4m-999m.json: $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/baseline_pluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Baseline/pluvial/P_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


$(temp)/flood_data/high_fluvial_1in1000_1m-2m.json: $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_fluvial_1in1000_1m-2m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_fluvial_1in1000_2m-3m.json: $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_fluvial_1in1000_2m-3m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_fluvial_1in1000_3m-4m.json: $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_fluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_fluvial_1in1000_4m-999m.json: $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_fluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/fluvial/FU_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


$(temp)/flood_data/high_pluvial_1in1000_1m-2m.json: $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_1m-2m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_pluvial_1in1000_1m-2m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_1m-2m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_pluvial_1in1000_2m-3m.json: $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_2m-3m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_pluvial_1in1000_2m-3m.json -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_2m-3m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_pluvial_1in1000_3m-4m.json: $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_3m-4m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_pluvial_1in1000_3m-4m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_3m-4m_threshold.shp -s_srs EPSG:4326

$(temp)/flood_data/high_pluvial_1in1000_4m-999m.json: $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_4m-999m_threshold.shp
	ogr2ogr -f GeoJSON $(temp)/flood_data/high_pluvial_1in1000_4m-999m.json  -t_srs EPSG:4326 $(in)/flood_data/FATHOM/Future_High/pluvial/P_1in1000_4m-999m_threshold.shp -s_srs EPSG:4326


$(out)/flood.mbtiles: $(temp)/flood_data/baseline_fluvial_1in1000_1m-2m.json $(temp)/flood_data/baseline_fluvial_1in1000_2m-3m.json $(temp)/flood_data/baseline_fluvial_1in1000_3m-4m.json $(temp)/flood_data/baseline_fluvial_1in1000_4m-999m.json $(temp)/flood_data/baseline_pluvial_1in1000_1m-2m.json $(temp)/flood_data/baseline_pluvial_1in1000_2m-3m.json $(temp)/flood_data/baseline_pluvial_1in1000_3m-4m.json $(temp)/flood_data/baseline_pluvial_1in1000_4m-999m.json $(temp)/flood_data/high_fluvial_1in1000_1m-2m.json $(temp)/flood_data/high_fluvial_1in1000_2m-3m.json $(temp)/flood_data/high_fluvial_1in1000_3m-4m.json $(temp)/flood_data/high_fluvial_1in1000_4m-999m.json $(temp)/flood_data/high_pluvial_1in1000_1m-2m.json $(temp)/flood_data/high_pluvial_1in1000_2m-3m.json $(temp)/flood_data/high_pluvial_1in1000_3m-4m.json $(temp)/flood_data/high_pluvial_1in1000_4m-999m.json
	rm -f $(out)/flood.mbtiles
	tippecanoe \
		-zg \
		--no-feature-limit \
		--no-line-simplification \
		--no-tile-size-limit \
		--extend-zooms-if-still-dropping \
		-o $(out)/flood.mbtiles \
		$(temp)/flood_data/*

clean:
	rm -f $(out)/*.mbtiles
