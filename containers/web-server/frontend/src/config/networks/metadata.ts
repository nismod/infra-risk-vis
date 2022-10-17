import { makeConfig } from '@/lib/helpers';

import { COLORS } from '@/config/colors';

/* Line widths:

-) elec_edges_high,
base: 1,
stops: [
  [7, 1],
  [12, 2],
  [16, 6],
],

-) elec_edges_low, pot_edges
base: 0.5,
stops: [
  [7, 0.5],
  [12, 1],
  [16, 3],
],

-) rail_edges
base: 1.5,
stops: [
  [7, 1.5],
  [12, 2],
  [16, 6],
],

-) road_edges
base: 0.5,
stops: [
  [7, 1.5],
  [12, 2],
  [16, 6],
],

-) all circle layers
base: 1.5,
stops: [
  [7, 3],
  [12, 4],
  [16, 12],
],

*/

export const NETWORKS_METADATA = makeConfig([
  {
    id: 'elec_edges_high',
    type: 'line',
    label: 'Power Lines (High Voltage)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_edges_low',
    type: 'line',
    label: 'Power Lines (Low Voltage)',
    color: COLORS.electricity_low.css,
  },
  {
    id: 'elec_nodes_diesel',
    type: 'circle',
    label: 'Power Generation (Diesel)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_gas',
    type: 'circle',
    label: 'Power Generation (Gas)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_hydro',
    type: 'circle',
    label: 'Power Generation (Hydro)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_solar',
    type: 'circle',
    label: 'Power Generation (Solar)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_wind',
    type: 'circle',
    label: 'Power Generation (Wind)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_pole',
    type: 'circle',
    label: 'Power Transmission (Poles)',
    color: COLORS.electricity_unknown.css,
  },
  {
    id: 'elec_nodes_substation',
    type: 'circle',
    label: 'Power Transmission (Substations)',
    color: COLORS.electricity_unknown.css,
  },
  {
    id: 'elec_nodes_demand',
    type: 'circle',
    label: 'Power Demand',
    color: COLORS.electricity_low.css,
  },
  {
    id: 'rail_edges',
    type: 'line',
    label: 'Railways',
    color: COLORS.railway.css,
  },
  {
    id: 'rail_stations',
    type: 'circle',
    label: 'Stations',
    color: COLORS.railway.css,
  },
  {
    id: 'rail_junctions',
    type: 'circle',
    label: 'Railway Junctions',
    color: COLORS.railway.css,
  },
  {
    id: 'road_edges_class_a',
    type: 'line',
    label: 'Roads (Class A)',
    color: COLORS.roads_class_a.css,
  },
  {
    id: 'road_edges_class_b',
    type: 'line',
    label: 'Roads (Class B)',
    color: COLORS.roads_class_b.css,
  },
  {
    id: 'road_edges_class_c',
    type: 'line',
    label: 'Roads (Class C)',
    color: COLORS.roads_class_c.css,
  },
  {
    id: 'road_edges_metro',
    type: 'line',
    label: 'Roads (Metro)',
    color: COLORS.roads_metro.css,
  },
  {
    id: 'road_edges_track',
    type: 'line',
    label: 'Roads (Tracks)',
    color: COLORS.roads_unknown.css,
  },
  {
    id: 'road_edges_other',
    type: 'line',
    label: 'Roads (Other)',
    color: COLORS.roads_unknown.css,
  },
  {
    id: 'road_bridges',
    type: 'circle',
    label: 'Bridges',
    color: COLORS.bridges.css,
  },
  {
    id: 'airport_runways',
    type: 'polygon',
    label: 'Airports (Runway)',
    color: COLORS.airports.css,
  },
  {
    id: 'airport_terminals',
    type: 'polygon',
    label: 'Airports (Terminal)',
    color: COLORS.airports.css,
  },
  {
    id: 'port_areas_break',
    type: 'polygon',
    label: 'Ports (Break)',
    color: COLORS.ports.css,
  },
  {
    id: 'port_areas_container',
    type: 'polygon',
    label: 'Ports (Container)',
    color: COLORS.ports.css,
  },
  {
    id: 'port_areas_industry',
    type: 'polygon',
    label: 'Ports (Industry)',
    color: COLORS.ports.css,
  },
  {
    id: 'port_areas_silo',
    type: 'polygon',
    label: 'Ports (Silo)',
    color: COLORS.ports.css,
  },
  {
    id: 'water_potable_edges',
    type: 'line',
    label: 'Water Supply Pipelines',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_booster',
    label: 'Water Supply (Booster Station)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_catchment',
    label: 'Water Supply (Catchment)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_entombment',
    label: 'Water Supply (Entombment)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_filter',
    label: 'Water Supply (Filter Plant)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_intake',
    label: 'Water Supply (Intake)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_well',
    label: 'Water Supply (Production Well)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_pump',
    label: 'Water Supply (Pump Station)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_relift',
    label: 'Water Supply (Relift Station)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_reservoir',
    label: 'Water Supply (Reservoir)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_river_source',
    label: 'Water Supply (River Source)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_spring',
    label: 'Water Supply (Spring)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_tank',
    label: 'Water Supply (Storage Tank)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_sump',
    label: 'Water Supply (Sump)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_tp',
    label: 'Water Supply (Treatment Plant)',
    type: 'circle',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_irrigation_edges',
    type: 'line',
    label: 'Irrigation Canals',
    color: COLORS.water_irrigation.css,
  },
  {
    id: 'water_irrigation_nodes',
    type: 'circle',
    label: 'Irrigation facilities',
    color: COLORS.water_irrigation.css,
  },
  {
    id: 'water_waste_sewer_gravity',
    type: 'line',
    label: 'Wastewater Pipelines (Gravity)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_sewer_pressure',
    type: 'line',
    label: 'Wastewater Pipelines (Pressure)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_sump',
    type: 'circle',
    label: 'Wastewater (Sump)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_pump',
    type: 'circle',
    label: 'Wastewater (Pump Station)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_relift',
    type: 'circle',
    label: 'Wastewater (Relift Station)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_wwtp',
    type: 'circle',
    label: 'Wastewater (Treament Plant)',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'buildings_commercial',
    type: 'polygon',
    label: 'Buildings (Commercial)',
    color: COLORS.buildings_commercial.css,
  },
  {
    id: 'buildings_industrial',
    type: 'polygon',
    label: 'Buildings (Industrial)',
    color: COLORS.buildings_industrial.css,
  },
  {
    id: 'buildings_institutional',
    type: 'polygon',
    label: 'Buildings (Institutional)',
    color: COLORS.buildings_institutional.css,
  },
  {
    id: 'buildings_mixed',
    type: 'polygon',
    label: 'Buildings (Mixed Use)',
    color: COLORS.buildings_mixed.css,
  },
  {
    id: 'buildings_other',
    type: 'polygon',
    label: 'Buildings (Other)',
    color: COLORS.buildings_other.css,
  },
  {
    id: 'buildings_recreation',
    type: 'polygon',
    label: 'Buildings (Recreation)',
    color: COLORS.buildings_recreation.css,
  },
  {
    id: 'buildings_residential',
    type: 'polygon',
    label: 'Buildings (Residential)',
    color: COLORS.buildings_residential.css,
  },
  {
    id: 'buildings_resort',
    type: 'polygon',
    label: 'Buildings (Resort)',
    color: COLORS.buildings_resort.css,
  },
]);
