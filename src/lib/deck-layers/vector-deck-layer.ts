import { MVTLayer, GeoJsonLayer } from 'deck.gl';

import { tileSelectionLayer, TileSelectionLayerOptions } from 'lib/deck-layers/tile-selection-layer';
import { mergeUpdateTriggers } from 'lib/deck-layers/utils';

export function vectorDeckLayer(selectionLayerParams: TileSelectionLayerOptions, ...props) {
  return new MVTLayer(
    {
      binary: false,
      autoHighlight: true,
      highlightColor: [0, 255, 255, 255],
      refinementStrategy: 'best-available',
      renderSubLayers: (tileProps) => [
        new GeoJsonLayer(tileProps),
        tileSelectionLayer(tileProps, selectionLayerParams),
      ],
    } as any,
    ...props,
    mergeUpdateTriggers(...props, {
      updateTriggers: {
        renderSubLayers: [selectionLayerParams.selectedFeatureId],
      },
    }),
  );
}
