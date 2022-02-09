import { LayerTree } from 'lib/layer-tree';
import { InteractionTarget } from 'lib/map/interactions/use-interactions';
import { ViewLayer, viewOnlyLayer } from 'lib/view-layers';

import { ViewLayerParams } from 'map/get-view-layers-spec';

import { BackgroundName } from './backgrounds';
import { boundaryLabelsLayer, BoundaryLevel } from './deck-layers/boundaries-layer';
import { labelsLayer } from './deck-layers/labels-layer';
import { VIEW_LAYERS } from './view-layers';
import { boundariesViewLayer } from './view-layers/boundaries-view-layer';

function getDataViewLayers(dataLayersSpec: Record<string, any>) {
  return Object.entries(dataLayersSpec).map(([viewLayerName, viewLayerParams]) => {
    const viewLayerConfig = VIEW_LAYERS[viewLayerName];
    const anyVisible = Object.values(viewLayerParams.visibility).some((x) => x);

    return anyVisible ? viewLayerConfig : null;
  });
}

export function getViewLayers(
  viewLayersSpec: Record<string, ViewLayerParams>,
  styleParams: any,
  selectedAsset: InteractionTarget<any>,
  showLabels: boolean,
  showBoundaries: boolean,
  boundaryLevel: BoundaryLevel,
  isRetina: boolean,
  background: BackgroundName,
): LayerTree<ViewLayer> {
  return [
    // boundaries layer
    showBoundaries && boundariesViewLayer(boundaryLevel),

    // all data layers
    getDataViewLayers(viewLayersSpec),

    showLabels && [
      // basemap labels layer
      viewOnlyLayer('labels', () => labelsLayer(isRetina)),

      // boundary labels layer
      showBoundaries &&
        viewOnlyLayer(`boundaries_${boundaryLevel}-text`, () => boundaryLabelsLayer(boundaryLevel, background)),
    ],
  ];
}
