import { MVTLayer, GeoJsonLayer } from 'deck.gl';

import { tileSelectionLayer, TileSelectionLayerOptions } from 'lib/deck-layers/tile-selection-layer';
import { mergeDeckProps } from './merge-props';

export function vectorDeckLayer(selectionLayerParams: TileSelectionLayerOptions, ...props) {
  return new MVTLayer(
    mergeDeckProps(
      {
        binary: false,
        autoHighlight: true,
        highlightColor: [0, 255, 255, 255],
        refinementStrategy: 'best-available',
        renderSubLayers: (tileProps) => [
          new GeoJsonLayer(tileProps),
          tileSelectionLayer(tileProps, selectionLayerParams),
        ],
        updateTriggers: {
          renderSubLayers: [selectionLayerParams.selectedFeatureId],
        },
      },
      props,
    ),
  );
}
