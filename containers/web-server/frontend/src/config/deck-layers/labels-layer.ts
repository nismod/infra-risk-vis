import GL from '@luma.gl/constants';

import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

export function labelsLayer(isRetina: boolean) {
  const scale = isRetina ? '@2x' : '';

  return rasterTileLayer(
    {
      [GL.TEXTURE_MAG_FILTER]: GL.NEAREST,
      transparentColor: [255, 255, 255, 0],
    },
    {
      id: 'labels',
      tileSize: 256,
      data: [
        `https://a.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}${scale}.png`,
        `https://b.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}${scale}.png`,
        `https://c.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}${scale}.png`,
        `https://d.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}${scale}.png`,
      ],
      refinementStrategy: 'no-overlap',
      getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 - 2000],
    },
  );
}
