import { InteractionTarget } from 'lib/map/interactions/use-interactions';

import { BackgroundName } from './backgrounds';
import { boundariesLayer, boundaryLabelsLayer, BoundaryLevel } from './deck-layers/boundaries-layer';
import { labelsLayer } from './deck-layers/labels-layer';
import { VIEW_LAYERS } from './view-layers';

function getDataDeckLayers(
  viewLayersSpec: Record<string, any>,
  zoom: number,
  styleParams: any,
  selectedFeature: InteractionTarget<any>,
) {
  return Object.entries(viewLayersSpec).map(([viewLayerName, viewLayerParams]) => {
    const viewLayerConfig = VIEW_LAYERS[viewLayerName];
    const anyVisible = Object.values(viewLayerParams.visibility).some((x) => x);

    let props: any = {
      id: viewLayerName,
      pickable: true,
      visible: anyVisible,
      minZoom: 3,
      maxZoom: 20,
    };

    let selectedFeatureId;
    if (selectedFeature?.viewLayer === viewLayerName) {
      selectedFeatureId = selectedFeature.target.id;
    }

    if (anyVisible) {
      return viewLayerConfig.fn({ props, ...viewLayerParams, zoom, styleParams, selectedFeatureId });
    }

    return null;
  });
}

export function getDeckLayers(
  viewLayersSpec: Record<string, any>,
  zoom: number,
  styleParams: any,
  selectedFeature: InteractionTarget<any>,
  showLabels: boolean,
  showBoundaries: boolean,
  boundaryLevel: BoundaryLevel,
  isRetina: boolean,
  background: BackgroundName,
) {
  return [
    // boundaries layer
    showBoundaries && boundariesLayer(boundaryLevel),

    // all data layers
    getDataDeckLayers(viewLayersSpec, zoom, styleParams, selectedFeature),

    showLabels && [
      // basemap labels layer
      labelsLayer(isRetina),

      // boundary labels layer
      showBoundaries && boundaryLabelsLayer(boundaryLevel, background),
    ],
  ];
}
