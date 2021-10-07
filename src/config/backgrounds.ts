import { makeConfig } from '../helpers';

export const backgroundConfig = makeConfig([
  {
    id: 'light',
    label: 'Map',
    source: {
      id: 'light',
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
    layer: {
      id: 'bg-light',
      source: 'light',
      type: 'raster',
    },
  },
  {
    id: 'satellite',
    label: 'Satellite',
    source: {
      id: 'satellite',
      type: 'raster',
      url: 'mapbox://mapbox.satellite',
      tileSize: 256,
    },
    layer: {
      id: 'bg-satellite',
      source: 'satellite',
      type: 'raster',
      'source-layer': 'mapbox_satellite_full',
    },
  },
]);

export type BackgroundName = keyof typeof backgroundConfig;
