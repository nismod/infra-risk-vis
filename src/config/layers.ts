export interface LayerDefinition {
  label: string;
  linear: boolean;
  color: string;
  style: object;
}

export const layers = {
  elec_edges_high: {
    linear: true,
    label: 'Power Lines (High Voltage)',
    color: '#eca926',
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
    linear: true,
    label: 'Power Lines (Low Voltage)',
    color: '#f1d75c',
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
    linear: false,
    label: 'Power Nodes',
    color: '#eca926',
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
    linear: true,
    label: 'Railways',
    color: '#444',
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
    linear: false,
    label: 'Stations',
    color: '#444',
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
    linear: true,
    label: 'Roads',
    color: '#b2afa',
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
    linear: false,
    label: 'Bridges',
    color: '#487dbc',
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
    linear: true,
    label: 'Water Supply Network',
    color: '#314386',
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
    linear: false,
    label: 'Water Abstraction',
    color: '#4d49bc',
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
  flood_fluvial_20: {
    linear: false,
    label: 'Fluvial flood (20 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_20',
      type: 'raster',
      source: 'flood_fluvial_20',
    },
  },
  flood_fluvial_50: {
    linear: false,
    label: 'Fluvial flood (50 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_50',
      type: 'raster',
      source: 'flood_fluvial_50',
    },
  },
  flood_fluvial_100: {
    linear: false,
    label: 'Fluvial flood (100 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_100',
      type: 'raster',
      source: 'flood_fluvial_100',
    },
  },
  flood_fluvial_200: {
    linear: false,
    label: 'Fluvial flood (200 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_200',
      type: 'raster',
      source: 'flood_fluvial_200',
    },
  },
  flood_fluvial_500: {
    linear: false,
    label: 'Fluvial flood (500 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_500',
      type: 'raster',
      source: 'flood_fluvial_500',
    },
  },
  flood_fluvial_1500: {
    linear: false,
    label: 'Fluvial flood (1500 year RP)',
    color: '#aaaaaa',
    style: {
      id: 'flood_fluvial_1500',
      type: 'raster',
      source: 'flood_fluvial_1500',
    },
  },
};

export type LayerName = keyof typeof layers;
