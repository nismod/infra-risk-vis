import { geoJsonLayer, mvtLayer } from './base';
import { DataLoaderOptions, dataLoaderLayer } from './data-loader-layer';
import { TileSelectionLayerOptions, tileSelectionLayer } from './tile-selection-layer';

interface SelectableMvtLayerOptions {
  selectionOptions: TileSelectionLayerOptions;
  dataLoaderOptions?: DataLoaderOptions;
}

export function selectableMvtLayer(
  { selectionOptions, dataLoaderOptions }: SelectableMvtLayerOptions,
  ...props
) {
  return mvtLayer(
    {
      binary: false,
      autoHighlight: true,
      highlightColor: [0, 255, 255, 255],
      refinementStrategy: 'best-available',
      renderSubLayers: (tileProps) => [
        geoJsonLayer(tileProps),
        tileSelectionLayer(tileProps, selectionOptions),
        dataLoaderOptions && dataLoaderLayer(tileProps, dataLoaderOptions),
      ],
      updateTriggers: {
        renderSubLayers: [selectionOptions.selectedFeatureId],
      },
    },
    props,
  );
}
