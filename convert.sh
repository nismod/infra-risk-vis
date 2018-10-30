ogr2ogr -f GeoJSON intermediate_data/admin_0_boundaries.json -t_srs EPSG:4326 incoming_data/boundaries/admin_0_boundaries.shp
ogr2ogr -f GeoJSON intermediate_data/admin_1_boundaries.json -t_srs EPSG:4326 incoming_data/boundaries/admin_1_boundaries.shp
ogr2ogr -f GeoJSON intermediate_data/economic_od_centroids.json -t_srs EPSG:4326 incoming_data/boundaries/economic_od_centroids.shp
ogr2ogr -f GeoJSON intermediate_data/economic_od_zones.json -t_srs EPSG:4326 incoming_data/boundaries/economic_od_zones.shp
ogr2ogr -f GeoJSON intermediate_data/physical_lakes.json -t_srs EPSG:4326 incoming_data/boundaries/physical_lakes.shp

ogr2ogr -f GeoJSON intermediate_data/air_edges.json -t_srs EPSG:4326 incoming_data/network/air_edges.shp
ogr2ogr -f GeoJSON intermediate_data/rail_edges.json -t_srs EPSG:4326 incoming_data/network/rail_edges.shp
ogr2ogr -f GeoJSON intermediate_data/road_edges_national.json -t_srs EPSG:4326 incoming_data/network/road_edges_national.shp
ogr2ogr -f GeoJSON intermediate_data/road_edges_provincial.json -t_srs EPSG:4326 incoming_data/network/road_edges_provincial.shp
ogr2ogr -f GeoJSON intermediate_data/road_edges_rural.json -t_srs EPSG:4326 incoming_data/network/road_edges_rural.shp
ogr2ogr -f GeoJSON intermediate_data/water_edges.json -t_srs EPSG:4326 incoming_data/network/water_edges.shp
ogr2ogr -f GeoJSON intermediate_data/water_nodes.json -t_srs EPSG:4326 incoming_data/network/water_nodes.shp

tippecanoe -o data/argentinia.mbtiles intermediate_data/*