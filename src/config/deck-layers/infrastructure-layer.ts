import { MVTLayer, GeoJsonLayer } from 'deck.gl';
import { tileSelectionLayer } from './tile-selection-layer';
import { mergeUpdateTriggers } from './utils';

export function infrastructureLayer({ selectedFeatureId }, ...props) {
  return new MVTLayer(
    {
      binary: false,
      autoHighlight: true,
      highlightColor: [0, 255, 255, 255],
      refinementStrategy: 'best-available',
      renderSubLayers: (tileProps) => [new GeoJsonLayer(tileProps), tileSelectionLayer(tileProps, selectedFeatureId)],
    } as any,
    ...props,
    mergeUpdateTriggers(...props, {
      updateTriggers: {
        renderSubLayers: [selectedFeatureId],
      },
    }),
  );
}
