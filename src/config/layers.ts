import { titleCase } from 'vega-lite';
import { makeConfig } from '../helpers';
import { COLORS } from './colors';

export interface LayerDefinition {
  deckLayer: string | { baseName: string; params: any };
  deckLayerParams?: any;
  label: string;
  type: string; //'line' | 'circle' | 'raster';
  color: string;
  getId?: (x) => string;
}

function floodLayer<F extends 'fluvial' | 'coastal', R extends number>(
  floodType: F,
  returnPeriod: R,
): LayerDefinition & { id: `hazard_${typeof floodType}_${typeof returnPeriod}` } {
  const id = `hazard_${floodType}_${returnPeriod}` as const;
  return {
    id,
    deckLayer: { baseName: 'hazard', params: { floodType, returnPeriod } },
    type: 'raster',
    label: `${titleCase(floodType)} flooding (${returnPeriod} year RP)`,
    color: '#aaaaaa',
    getId: ({ floodType, returnPeriod }) => `hazard_${floodType}_${returnPeriod}`,
  };
}

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
export const LAYERS = makeConfig([
  {
    id: 'elec_edges_high',
    deckLayer: 'elec_edges',
    type: 'line',
    label: 'Power Lines (High Voltage)',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'elec_edges_low',
    deckLayer: 'elec_edges',
    type: 'line',
    label: 'Power Lines (Low Voltage)',
    color: COLORS.electricity_low.css,
  },
  {
    id: 'elec_nodes',
    deckLayer: 'elec_nodes',
    type: 'circle',
    label: 'Power Nodes',
    color: COLORS.electricity_high.css,
  },
  {
    id: 'rail_edges',
    deckLayer: 'rail_edges',
    type: 'line',
    label: 'Railways',
    color: COLORS.railway.css,
  },
  {
    id: 'rail_nodes',
    deckLayer: 'rail_nodes',
    type: 'circle',
    label: 'Stations',
    color: COLORS.railway.css,
  },
  {
    id: 'road_edges',
    deckLayer: 'road_edges',
    type: 'line',
    label: 'Roads',
    color: COLORS.roads_unknown.css,
  },
  {
    id: 'bridges',
    deckLayer: 'bridges',
    type: 'circle',
    label: 'Bridges',
    color: COLORS.bridges.css,
  },
  {
    id: 'pot_edges',
    deckLayer: 'pot_edges',
    type: 'line',
    label: 'Water Supply Network',
    color: COLORS.water_edges.css,
  },
  {
    id: 'abs_nodes',
    deckLayer: 'abs_nodes',
    type: 'circle',
    label: 'Water Abstraction',
    color: COLORS.water_abstraction.css,
  },
  floodLayer('fluvial', 20),
  floodLayer('fluvial', 50),
  floodLayer('fluvial', 100),
  floodLayer('fluvial', 200),
  floodLayer('fluvial', 500),
  floodLayer('fluvial', 1500),
  floodLayer('coastal', 1),
  floodLayer('coastal', 2),
  floodLayer('coastal', 5),
  floodLayer('coastal', 10),
  floodLayer('coastal', 50),
  floodLayer('coastal', 100),
]);

export type LayerName = keyof typeof LAYERS;
