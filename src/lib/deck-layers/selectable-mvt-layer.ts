import { tileSelectionLayer, TileSelectionLayerOptions } from './tile-selection-layer';
import { geoJsonLayer, mvtLayer } from './base';

export function selectableMvtLayer(selectionLayerParams: TileSelectionLayerOptions, ...props) {
  return mvtLayer(
    {
      binary: false,
      autoHighlight: true,
      highlightColor: [0, 255, 255, 255],
      refinementStrategy: 'best-available',
      renderSubLayers: (tileProps) => [geoJsonLayer(tileProps), tileSelectionLayer(tileProps, selectionLayerParams)],
      updateTriggers: {
        renderSubLayers: [selectionLayerParams.selectedFeatureId],
      },
    },
    props,
  );
}
