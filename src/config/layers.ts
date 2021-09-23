import { titleCase } from 'vega-lite';

export interface LayerDefinition {
  label: string;
  type: 'line' | 'circle' | 'raster';
  color: string;
  style: object;
  sourceUrl: string;
  width?: any;
  radius?: any;
}

const rasterColormaps = {
  fluvial: 'blues',
  coastal: 'greens',
};

const rasterColormapRanges = {
  fluvial: '[0,10]',
  coastal: '[0,3.5]',
};

function floodLayer(floodType: string, returnPeriod: number): LayerDefinition {
  const id = `flood_${floodType}_${returnPeriod}`;
  return {
    type: 'raster',
    label: `${titleCase(floodType)} flooding (${returnPeriod} year RP)`,
    color: '#aaaaaa',
    style: {
      id,
      type: 'raster',
      source: id,
      paint: {
        'raster-fade-duration': 0,
      },
    },
    sourceUrl: `http://localhost:5000/singleband/${floodType}/${returnPeriod}/raw/{z}/{x}/{y}.png?colormap=${rasterColormaps[floodType]}&stretch_range=${rasterColormapRanges[floodType]}`,
  };
}

const LINE_MIN_PIXELS = 1;
const POINT_MIN_PIXELS = 3;

const layersConfig = {
  elec_edges_high: {
    type: 'line',
    label: 'Power Lines (High Voltage)',
    color: '#eca926',
    sourceUrl: 'http://localhost:8080/data/elec_edges.json',
    width: {
      value: 5,
      unit: 'meters',
      minPixels: LINE_MIN_PIXELS,
      maxPixels: 10,
    },
    style: {
      id: 'elec_edges_high',
      type: 'line',
      source: 'elec_edges',
      'source-layer': 'elec_edges',
      minzoom: 3,
      filter: ['==', 'High Voltage', ['get', 'line_type']],
      layout: {
        'line-cap': 'round',
        'line-join': 'round',
      },
      paint: {
        'line-color': '#eca926',
        'line-width': {
          base: 1,
          stops: [
            [7, 1],
            [12, 2],
            [16, 6],
          ],
        },
      },
    },
  },
  elec_edges_low: {
    type: 'line',
    label: 'Power Lines (Low Voltage)',
    color: '#f1d75c',
    sourceUrl: 'http://localhost:8080/data/elec_edges.json',
    width: {
      value: 5,
      unit: 'meters',
      minPixels: LINE_MIN_PIXELS,
      maxPixels: 10,
    },
    style: {
      id: 'elec_edges_low',
      type: 'line',
      source: 'elec_edges',
      'source-layer': 'elec_edges',
      minzoom: 3,
      filter: ['==', 'Low Voltage', ['get', 'line_type']],
      layout: {
        'line-cap': 'round',
        'line-join': 'round',
      },
      paint: {
        'line-color': '#f1d75c',
        'line-width': {
          base: 0.5,
          stops: [
            [7, 0.5],
            [12, 1],
            [16, 3],
          ],
        },
      },
    },
  },
  elec_nodes: {
    type: 'circle',
    label: 'Power Nodes',
    color: '#eca926',
    sourceUrl: 'http://localhost:8080/data/elec_nodes.json',
    radius: {
      value: 10,
      unit: 'meters',
      minPixels: POINT_MIN_PIXELS,
      maxPixels: 20,
    },
    style: {
      id: 'elec_nodes',
      type: 'circle',
      source: 'elec_nodes',
      'source-layer': 'elec_nodes',
      minzoom: 3,
      paint: {
        'circle-color': '#eca926',
        'circle-radius': {
          base: 1.5,
          stops: [
            [7, 3],
            [12, 4],
            [16, 12],
          ],
        },
      },
    },
  },
  rail_edges: {
    type: 'line',
    label: 'Railways',
    color: '#444',
    sourceUrl: 'http://localhost:8080/data/rail_edges.json',
    width: {
      value: 5,
      unit: 'meters',
      minPixels: LINE_MIN_PIXELS,
      maxPixels: 10,
    },
    style: {
      id: 'rail_edges',
      type: 'line',
      source: 'rail_edges',
      'source-layer': 'rail_edges',
      minzoom: 3,
      layout: {
        'line-cap': 'round',
        'line-join': 'round',
      },
      paint: {
        'line-color': '#444',
        'line-width': {
          base: 1.5,
          stops: [
            [7, 1.5],
            [12, 2],
            [16, 6],
          ],
        },
      },
    },
  },
  rail_nodes: {
    type: 'circle',
    label: 'Stations',
    color: '#444',
    sourceUrl: 'http://localhost:8080/data/rail_nodes.json',
    radius: {
      value: 10,
      unit: 'meters',
      minPixels: POINT_MIN_PIXELS,
      maxPixels: 20,
    },
    style: {
      id: 'rail_nodes',
      type: 'circle',
      source: 'rail_nodes',
      'source-layer': 'rail_nodes',
      minzoom: 3,
      paint: {
        'circle-color': '#444',
        'circle-radius': {
          base: 1.5,
          stops: [
            [7, 3],
            [12, 4],
            [16, 12],
          ],
        },
      },
    },
  },
  road_edges: {
    type: 'line',
    label: 'Roads',
    color: '#b2afaa',
    sourceUrl: 'http://localhost:8080/data/road_edges.json',
    width: {
      value: 5,
      unit: 'meters',
      minPixels: LINE_MIN_PIXELS,
      maxPixels: 10,
    },
    style: {
      id: 'road_edges',
      type: 'line',
      source: 'road_edges',
      'source-layer': 'road_edges',
      minzoom: 3,
      layout: {
        'line-cap': 'round',
        'line-join': 'round',
      },
      paint: {
        'line-color': [
          'match',
          ['get', 'class'],
          'CLASS A',
          '#941339',
          'CLASS B',
          '#cb3e4e',
          'CLASS C',
          '#8471a8',
          'METRO',
          '#487dbc',
          '#b2afaa',
        ],
        'line-width': {
          base: 0.5,
          stops: [
            [7, 1.5],
            [12, 2],
            [16, 6],
          ],
        },
      },
    },
  },
  bridges: {
    type: 'circle',
    label: 'Bridges',
    color: '#487dbc',
    sourceUrl: 'http://localhost:8080/data/bridges.json',
    radius: {
      value: 10,
      unit: 'meters',
      minPixels: POINT_MIN_PIXELS,
      maxPixels: 20,
    },
    style: {
      id: 'bridges',
      type: 'circle',
      source: 'bridges',
      'source-layer': 'bridges',
      minzoom: 3,
      paint: {
        'circle-color': '#487dbc',
        'circle-radius': {
          base: 1.5,
          stops: [
            [7, 3],
            [12, 4],
            [16, 12],
          ],
        },
      },
    },
  },
  pot_edges: {
    type: 'line',
    label: 'Water Supply Network',
    color: '#314386',
    sourceUrl: 'http://localhost:8080/data/pot_edges.json',
    width: {
      value: 5,
      unit: 'meters',
      minPixels: LINE_MIN_PIXELS,
      maxPixels: 10,
    },
    style: {
      id: 'pot_edges',
      type: 'line',
      source: 'pot_edges',
      'source-layer': 'pot_edges',
      minzoom: 3,
      layout: {
        'line-cap': 'round',
        'line-join': 'round',
      },
      paint: {
        'line-color': '#314386',
        'line-width': {
          base: 0.5,
          stops: [
            [7, 0.5],
            [12, 1],
            [16, 3],
          ],
        },
      },
    },
  },
  abs_nodes: {
    type: 'circle',
    label: 'Water Abstraction',
    color: '#4d49bc',
    sourceUrl: 'http://localhost:8080/data/abs_nodes.json',
    radius: {
      value: 10,
      unit: 'meters',
      minPixels: POINT_MIN_PIXELS,
      maxPixels: 20,
    },
    style: {
      id: 'abs_nodes',
      type: 'circle',
      source: 'abs_nodes',
      'source-layer': 'abs_nodes',
      minzoom: 3,
      paint: {
        'circle-color': '#4d49bc',
        'circle-radius': {
          base: 1.5,
          stops: [
            [7, 3],
            [12, 4],
            [16, 12],
          ],
        },
      },
    },
  },
  flood_fluvial_20: floodLayer('fluvial', 20),
  flood_fluvial_50: floodLayer('fluvial', 50),
  flood_fluvial_100: floodLayer('fluvial', 100),
  flood_fluvial_200: floodLayer('fluvial', 200),
  flood_fluvial_500: floodLayer('fluvial', 500),
  flood_fluvial_1500: floodLayer('fluvial', 1500),
  flood_coastal_1: floodLayer('coastal', 1),
  flood_coastal_2: floodLayer('coastal', 2),
  flood_coastal_5: floodLayer('coastal', 5),
  flood_coastal_10: floodLayer('coastal', 10),
  flood_coastal_50: floodLayer('coastal', 50),
  flood_coastal_100: floodLayer('coastal', 100),
};

export type LayerName = keyof typeof layersConfig;

export const layers = layersConfig as Record<LayerName, LayerDefinition>;
