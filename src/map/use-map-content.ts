import { useMemo } from 'react';

/**
 * get map style and layers definition based on:
 * - selected background
 * - selected layers, filters
 * - selected data visualisation
 * - any highlights / selections
 */

interface MapParams {
  background: 'satellite' | 'light';
}

function getMapSources() {
  return {
    satellite: {
      type: 'raster',
      url: 'mapbox://mapbox.satellite',
      tileSize: 256,
    },
    light: {
      type: 'raster',
      tiles: [
        'https://tiles.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://a.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://b.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://c.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://d.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
      ],
      tileSize: 256,
    },

    road_edges: {
      type: 'vector',
      url: 'http://localhost:8080/data/road_edges.json',
    },
    bridges: {
      type: 'vector',
      url: 'http://localhost:8080/data/bridges.json',
    },
    elec_edges: {
      type: 'vector',
      url: 'http://localhost:8080/data/elec_edges.json',
    },
    elec_nodes: {
      type: 'vector',
      url: 'http://localhost:8080/data/elec_nodes.json',
    },
    rail_edges: {
      type: 'vector',
      url: 'http://localhost:8080/data/rail_edges.json',
    },
    rail_nodes: {
      type: 'vector',
      url: 'http://localhost:8080/data/rail_nodes.json',
    },
    pot_edges: {
      type: 'vector',
      url: 'http://localhost:8080/data/pot_edges.json',
    },
    abs_nodes: {
      type: 'vector',
      url: 'http://localhost:8080/data/abs_nodes.json',
    },
  };
}

const overviewLayers = [
  {
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
  {
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
  {
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
  {
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
  {
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
  {
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
  {
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
  {
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
  {
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
];

function getMapLayers({ background }: MapParams) {
  let res = [];

  res.push({
    id: 'bg-satellite',
    source: 'satellite',
    type: 'raster',
    'source-layer': 'mapbox_satellite_full',
    layout: {
      visibility: background === 'satellite' ? 'visible' : 'none',
    },
  });

  res.push({
    id: 'bg-light',
    source: 'light',
    type: 'raster',
    layout: {
      visibility: background === 'light' ? 'visible' : 'none',
    },
  });

  res = res.concat(overviewLayers);

  // for (let i = 0; i < res.length - 1; i++) {
  //   res[i].beforeId = res[i + 1].id;
  // }

  return res;
}

export function useMapContent(params: MapParams) {
  const sources = useMemo(() => getMapSources(), []);
  const layers = useMemo(() => getMapLayers(params), [params]);

  const res = useMemo(
    () => ({
      version: 8,
      sources,
      layers,
      glyphs: 'http://localhost:8080/{fontstack}/{range}.pbf',
    }),
    [sources, layers],
  );

  return res;
}
