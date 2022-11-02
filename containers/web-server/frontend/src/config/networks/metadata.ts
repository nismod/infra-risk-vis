import { makeConfig } from '@/lib/helpers';

import { INFRASTRUCTURE_COLORS } from '@/config/networks/colors';

import { AssetMetadata } from '../assets/metadata';

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

export const NETWORK_LAYERS = [
  'elec_edges_high',
  'elec_edges_low',
  'elec_nodes_diesel',
  'elec_nodes_gas',
  'elec_nodes_hydro',
  'elec_nodes_solar',
  'elec_nodes_wind',
  'elec_nodes_pole',
  'elec_nodes_substation',
  'elec_nodes_demand',
  'rail_edges',
  'rail_stations',
  'rail_nodes',
  'road_edges_class_a',
  'road_edges_class_b',
  'road_edges_class_c',
  'road_edges_metro',
  'road_edges_track',
  'road_edges_other',
  'road_bridges',
  'airport_runways',
  'airport_terminals',
  'port_areas_break',
  'port_areas_container',
  'port_areas_industry',
  'port_areas_silo',
  'water_potable_edges',
  'water_potable_nodes_booster',
  'water_potable_nodes_catchment',
  'water_potable_nodes_entombment',
  'water_potable_nodes_filter',
  'water_potable_nodes_intake',
  'water_potable_nodes_well',
  'water_potable_nodes_pump',
  'water_potable_nodes_relift',
  'water_potable_nodes_reservoir',
  'water_potable_nodes_river_source',
  'water_potable_nodes_spring',
  'water_potable_nodes_tank',
  'water_potable_nodes_sump',
  'water_potable_nodes_tp',
  'water_irrigation_edges',
  'water_irrigation_nodes',
  'water_waste_sewer_gravity',
  'water_waste_sewer_pressure',
  'water_waste_nodes_sump',
  'water_waste_nodes_pump',
  'water_waste_nodes_relift',
  'water_waste_nodes_wwtp',
] as const;

export type NetworkLayerType = typeof NETWORK_LAYERS[number];

export const NETWORKS_METADATA = makeConfig<AssetMetadata, NetworkLayerType>([
  {
    id: 'elec_edges_high',
    type: 'line',
    label: 'Power Lines (High Voltage)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_edges_low',
    type: 'line',
    label: 'Power Lines (Low Voltage)',
    color: INFRASTRUCTURE_COLORS.electricity_low.css,
  },
  {
    id: 'elec_nodes_diesel',
    type: 'circle',
    label: 'Power Generation (Diesel)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_gas',
    type: 'circle',
    label: 'Power Generation (Gas)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_hydro',
    type: 'circle',
    label: 'Power Generation (Hydro)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_solar',
    type: 'circle',
    label: 'Power Generation (Solar)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_wind',
    type: 'circle',
    label: 'Power Generation (Wind)',
    color: INFRASTRUCTURE_COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_pole',
    type: 'circle',
    label: 'Power Transmission (Poles)',
    color: INFRASTRUCTURE_COLORS.electricity_unknown.css,
  },
  {
    id: 'elec_nodes_substation',
    type: 'circle',
    label: 'Power Transmission (Substations)',
    color: INFRASTRUCTURE_COLORS.electricity_unknown.css,
  },
  {
    id: 'elec_nodes_demand',
    type: 'circle',
    label: 'Power Demand',
    color: INFRASTRUCTURE_COLORS.electricity_low.css,
  },
  {
    id: 'rail_edges',
    type: 'line',
    label: 'Railways',
    color: INFRASTRUCTURE_COLORS.railway.css,
  },
  {
    id: 'rail_stations',
    type: 'circle',
    label: 'Stations',
    color: INFRASTRUCTURE_COLORS.railway.css,
  },
  {
    id: 'rail_nodes',
    type: 'circle',
    label: 'Railway Junctions',
    color: INFRASTRUCTURE_COLORS.railway.css,
  },
  {
    id: 'road_edges_class_a',
    type: 'line',
    label: 'Roads (Motorway)',
    color: INFRASTRUCTURE_COLORS.roads_class_a.css,
  },
  {
    id: 'road_edges_class_b',
    type: 'line',
    label: 'Roads (Trunk)',
    color: INFRASTRUCTURE_COLORS.roads_class_b.css,
  },
  {
    id: 'road_edges_class_c',
    type: 'line',
    label: 'Roads (Primary)',
    color: INFRASTRUCTURE_COLORS.roads_class_c.css,
  },
  {
    id: 'road_edges_metro',
    type: 'line',
    label: 'Roads (Secondary)',
    color: INFRASTRUCTURE_COLORS.roads_metro.css,
  },
  {
    id: 'road_edges_track',
    type: 'line',
    label: 'Roads (Tertiary)',
    color: INFRASTRUCTURE_COLORS.roads_unknown.css,
  },
  {
    id: 'road_edges_other',
    type: 'line',
    label: 'Roads (Other)',
    color: INFRASTRUCTURE_COLORS.roads_unknown.css,
  },
  {
    id: 'road_bridges',
    type: 'circle',
    label: 'Bridges',
    color: INFRASTRUCTURE_COLORS.bridges.css,
  },
  {
    id: 'airport_runways',
    type: 'polygon',
    label: 'Airports (Runway)',
    color: INFRASTRUCTURE_COLORS.airports.css,
  },
  {
    id: 'airport_terminals',
    type: 'polygon',
    label: 'Airports (Terminal)',
    color: INFRASTRUCTURE_COLORS.airports.css,
  },
  {
    id: 'port_areas_break',
    type: 'polygon',
    label: 'Ports (Break)',
    color: INFRASTRUCTURE_COLORS.ports.css,
  },
  {
    id: 'port_areas_container',
    type: 'polygon',
    label: 'Ports (Container)',
    color: INFRASTRUCTURE_COLORS.ports.css,
  },
  {
    id: 'port_areas_industry',
    type: 'polygon',
    label: 'Ports (Industry)',
    color: INFRASTRUCTURE_COLORS.ports.css,
  },
  {
    id: 'port_areas_silo',
    type: 'polygon',
    label: 'Ports (Silo)',
    color: INFRASTRUCTURE_COLORS.ports.css,
  },
  {
    id: 'water_potable_edges',
    type: 'line',
    label: 'Water Supply Pipelines',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_booster',
    label: 'Water Supply (Booster Station)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_catchment',
    label: 'Water Supply (Catchment)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_entombment',
    label: 'Water Supply (Entombment)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_filter',
    label: 'Water Supply (Filter Plant)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_intake',
    label: 'Water Supply (Intake)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_well',
    label: 'Water Supply (Production Well)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_pump',
    label: 'Water Supply (Pump Station)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_relift',
    label: 'Water Supply (Relift Station)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_reservoir',
    label: 'Water Supply (Reservoir)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_river_source',
    label: 'Water Supply (River Source)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_spring',
    label: 'Water Supply (Spring)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_tank',
    label: 'Water Supply (Storage Tank)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_sump',
    label: 'Water Supply (Sump)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes_tp',
    label: 'Water Supply (Treatment Plant)',
    type: 'circle',
    color: INFRASTRUCTURE_COLORS.water_supply.css,
  },
  {
    id: 'water_irrigation_edges',
    type: 'line',
    label: 'Irrigation Canals',
    color: INFRASTRUCTURE_COLORS.water_irrigation.css,
  },
  {
    id: 'water_irrigation_nodes',
    type: 'circle',
    label: 'Irrigation facilities',
    color: INFRASTRUCTURE_COLORS.water_irrigation.css,
  },
  {
    id: 'water_waste_sewer_gravity',
    type: 'line',
    label: 'Wastewater Pipelines (Gravity)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_sewer_pressure',
    type: 'line',
    label: 'Wastewater Pipelines (Pressure)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_sump',
    type: 'circle',
    label: 'Wastewater (Sump)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_pump',
    type: 'circle',
    label: 'Wastewater (Pump Station)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_relift',
    type: 'circle',
    label: 'Wastewater (Relift Station)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes_wwtp',
    type: 'circle',
    label: 'Wastewater (Treament Plant)',
    color: INFRASTRUCTURE_COLORS.water_wastewater.css,
  },
]);
