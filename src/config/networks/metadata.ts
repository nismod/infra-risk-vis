import { COLORS } from 'config/colors';
import { makeConfig } from 'lib/helpers';

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
    id: 'elec_nodes_source',
    type: 'circle',
    label: 'Power Nodes (Generation)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_nodes_junction',
    type: 'circle',
    label: 'Power Nodes (Junctions)',
    color: COLORS.electricity_unknown.css,
  },
  {
    id: 'elec_nodes_sink',
    type: 'circle',
    label: 'Power Nodes (Demand)',
    color: COLORS.electricity_low.css,
  },
  {
    id: 'rail_edges',
    type: 'line',
    label: 'Railways',
    color: COLORS.railway.css,
  },
  {
    id: 'rail_nodes',
    type: 'circle',
    label: 'Stations',
    color: COLORS.railway.css,
  },
  {
    id: 'road_edges',
    type: 'line',
    label: 'Roads',
    color: COLORS.roads_unknown.css,
  },
  {
    id: 'road_bridges',
    type: 'circle',
    label: 'Bridges',
    color: COLORS.bridges.css,
  },
  {
    id: 'airport_areas',
    type: 'polygon',
    label: 'Airports',
    color: COLORS.airports.css,
  },
  {
    id: 'port_areas',
    type: 'polygon',
    label: 'Ports',
    color: COLORS.ports.css,
  },
  {
    id: 'water_potable_edges',
    type: 'line',
    label: 'Water Supply Pipelines',
    color: COLORS.water_supply.css,
  },
  {
    id: 'water_potable_nodes',
    type: 'circle',
    label: 'Water Supply Facilities',
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
    id: 'water_waste_edges',
    type: 'line',
    label: 'Wastewater Pipelines',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'water_waste_nodes',
    type: 'circle',
    label: 'Wastewater Facilities',
    color: COLORS.water_wastewater.css,
  },
  {
    id: 'buildings',
    type: 'square',
    label: 'Buildings',
    color: COLORS.buildings.css,
  },
]);

export const NETWORK_LAYER_NAMES = [
  'elec_edges_high',
  'elec_edges_low',
  'elec_nodes_source',
  'elec_nodes_sink',
  'elec_nodes_junction',
  'rail_edges',
  'rail_nodes',
  'road_edges',
  'road_bridges',
  'port_areas',
  'airport_areas',
  'water_potable_edges',
  'water_potable_nodes',
  'water_irrigation_edges',
  'water_irrigation_nodes',
  'water_waste_edges',
  'water_waste_nodes',
];
