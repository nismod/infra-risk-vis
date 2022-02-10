import { GeoJsonLayer } from 'deck.gl';
import { DataFilterExtension } from '@deck.gl/extensions';

import { mergeUpdateTriggers } from './utils';

export function tileSelectionLayer(tileProps, selectedFeatureId) {
  const layer = new GeoJsonLayer<{ id: any }>(
    tileProps,
    {
      id: tileProps.id + '-selection',
      pickable: false,
      getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 - 1000],
      visible: selectedFeatureId != null,

      getLineWidth: 2,
      lineWidthUnits: 'pixels',
      getFillColor: [0, 255, 255, 255],
      getLineColor: [0, 255, 255, 255],

      // use on-GPU filter extension to only show the selected feature
      getFilterValue: (x) => (x.id == selectedFeatureId ? 1 : 0),
      filterRange: [1, 1],
      extensions: [new DataFilterExtension({ filterSize: 1 })],
    } as any,
    mergeUpdateTriggers(tileProps, {
      updateTriggers: {
        getLineWidth: [selectedFeatureId],
        getFilterValue: [selectedFeatureId],
      },
    }),
  );
  return layer;
}
