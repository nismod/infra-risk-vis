import { useMemo } from 'react';

import { LayerName, layers } from '../config/layers';
import { ViewName, views } from '../config/views';

/**
 * get map style and layers definition based on:
 * - selected background
 * - selected layers, filters
 * - selected data visualisation
 * - any highlights / selections
 */

export interface MapParams {
  background: 'satellite' | 'light';
  view: ViewName;
  dataLayerSelection: Record<LayerName, boolean>;
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

function visible(isVisible: boolean): 'visible' | 'none' {
  return isVisible ? 'visible' : 'none';
}

function getMapLayers({ background, view, dataLayerSelection }: MapParams) {
  let res = [];

  res.push({
    id: 'bg-satellite',
    source: 'satellite',
    type: 'raster',
    'source-layer': 'mapbox_satellite_full',
    layout: {
      visibility: visible(background === 'satellite'),
    },
  });

  res.push({
    id: 'bg-light',
    source: 'light',
    type: 'raster',
    layout: {
      visibility: visible(background === 'light'),
    },
  });

  for (const layerName of views[view].layers) {
    const layerStyle = layers[layerName].style;

    // TODO: prevent overwriting layout from original style - perform deep merge instead
    res.push({ ...layerStyle, layout: { visibility: visible(dataLayerSelection[layerName]) } });
  }

  // TODO: handle setting beforeId if needed
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
