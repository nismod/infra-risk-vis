set -e
set -x
psql -c '\copy damages_expected from raw_data/processed_data/input/20221027_global_rail_stations_rail_nodes_exp_0.csv with csv header'
