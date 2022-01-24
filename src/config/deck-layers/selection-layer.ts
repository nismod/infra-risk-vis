import { GeoJsonLayer } from 'deck.gl';

import { lineStyle, pointRadius } from './utils';

export function selectionLayer(feature, zoom) {
  return new GeoJsonLayer<any>(
    {
      id: 'selection',
      data: [feature],
      getFillColor: [0, 255, 255],
      getLineColor: [0, 255, 255],
      pickable: false,
    },
    lineStyle(zoom),
    pointRadius(zoom),
  );
}
