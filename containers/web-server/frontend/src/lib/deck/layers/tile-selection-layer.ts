import { DataFilterExtension } from '@deck.gl/extensions/typed';

import { geoJsonLayer } from './base';

export interface TileSelectionLayerOptions {
  selectedFeatureId: number | null;
  selectionFillColor?: [number, number, number, number];
  selectionLineColor?: [number, number, number, number];
  polygonOffset?: number;
}
export function tileSelectionLayer(
  tileProps,
  {
    selectedFeatureId,
    selectionFillColor = [0, 255, 255, 255],
    selectionLineColor = [0, 255, 255, 255],
    polygonOffset = 0,
  }: TileSelectionLayerOptions,
) {
  return geoJsonLayer(tileProps, {
    id: tileProps.id + '-selection',
    pickable: false,
    getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 + polygonOffset],
    visible: selectedFeatureId != null,
    refinementStrategy: 'no-overlap',

    getLineWidth: 2,
    lineWidthUnits: 'pixels',
    getFillColor: selectionFillColor,
    getLineColor: selectionLineColor,

    updateTriggers: {
      getLineWidth: [selectedFeatureId],
      getFilterValue: [selectedFeatureId],
    },

    // use on-GPU filter extension to only show the selected feature
    getFilterValue: (x) => (x.id === selectedFeatureId ? 1 : 0),
    filterRange: [1, 1],
    extensions: [new DataFilterExtension({ filterSize: 1 })],
  });
}
