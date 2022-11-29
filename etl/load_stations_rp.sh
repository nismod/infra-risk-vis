set -e
set -x
psql -c '\copy damages_rp from raw_data/processed_data/input/20221027_global_rail_stations_rail_nodes_rp_0.csv with csv header'
