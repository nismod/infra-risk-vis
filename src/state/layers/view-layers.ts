import { boundaryLabelsLayer } from 'config/regions/boundaries-deck-layer';
import { boundariesViewLayer } from 'config/regions/boundaries-view-layer';
import { hazardViewLayer } from 'config/hazards/hazard-view-layer';
import { ViewLayer, viewOnlyLayer } from 'lib/data-map/view-layers';
import { backgroundState, boundaryLevelState, showBoundariesState, showLabelsState } from 'map/layers/layers-state';
import { selector } from 'recoil';
import { dataParamsByGroupState } from 'state/data-params';
import { hazardVisibilityState } from 'state/hazards/hazard-visibility';
import { networkSelectionState } from 'state/network-selection';
import { truthyKeys } from 'lib/helpers';
import { INFRASTRUCTURE_VIEW_LAYERS } from 'config/networks/view-layers';
import { labelsLayer } from 'config/deck-layers/labels-layer';
import { isRetinaState } from 'state/is-retina';
import { HazardParams } from 'config/hazards/domains';
import { LayerTree } from 'lib/layer-tree';

import { showPopulationState } from '../population';
import { populationViewLayer } from 'config/regions/population-view-layer';

const hazardLayerState = selector<ViewLayer[]>({
  key: 'hazardLayerState',
  get: ({ get }) =>
    truthyKeys(get(hazardVisibilityState)).map((hazard) =>
      hazardViewLayer(hazard, get(dataParamsByGroupState(hazard)) as HazardParams),
    ),
});

const networkLayersState = selector<ViewLayer[]>({
  key: 'networkLayersState',
  get: ({ get }) => truthyKeys(get(networkSelectionState)).map((network) => INFRASTRUCTURE_VIEW_LAYERS[network]),
});

export const viewLayersState = selector<LayerTree<ViewLayer>>({
  key: 'viewLayersState',
  get: ({ get }) => {
    const showBoundaries = get(showBoundariesState);
    const boundaryLevel = get(boundaryLevelState);
    const background = get(backgroundState);
    const showLabels = get(showLabelsState);
    const isRetina = get(isRetinaState);

    return [
      // administrative region boundaries or population density
      get(showPopulationState)
        ? populationViewLayer(boundaryLevel)
        : showBoundaries && boundariesViewLayer(boundaryLevel),

      // hazard data layers
      get(hazardLayerState),

      // network data layers
      get(networkLayersState),

      showLabels && [
        // basemap labels
        viewOnlyLayer('labels', () => labelsLayer(isRetina)),

        // administrative regions labels
        showBoundaries &&
          viewOnlyLayer(`boundaries_${boundaryLevel}-text`, () => boundaryLabelsLayer(boundaryLevel, background)),
      ],
    ];
  },
});
