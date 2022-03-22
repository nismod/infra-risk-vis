import { regionLabelsDeckLayer } from 'config/regions/region-labels-deck-layer';
import { regionBoundariesViewLayer } from 'config/regions/boundaries-view-layer';
import { hazardViewLayer } from 'config/hazards/hazard-view-layer';
import { ViewLayer, viewOnlyLayer } from 'lib/data-map/view-layers';
import { backgroundState, showLabelsState } from 'map/layers/layers-state';
import { selector } from 'recoil';
import { dataParamsByGroupState } from 'state/data-params';
import { hazardVisibilityState } from 'state/hazards/hazard-visibility';
import { networkSelectionState } from 'state/networks/network-selection';
import { truthyKeys } from 'lib/helpers';
import { INFRASTRUCTURE_VIEW_LAYERS } from 'config/networks/view-layers';
import { labelsLayer } from 'config/deck-layers/labels-layer';
import { isRetinaState } from 'state/is-retina';
import { HazardParams } from 'config/hazards/domains';
import { ConfigTree } from 'lib/config-tree';

import { populationViewLayer } from 'config/regions/population-view-layer';
import { regionLevelState, showPopulationState } from 'state/regions';
import { sectionVisibilityState } from 'state/sections';
import { buildingsViewLayer } from 'config/buildings/buildings-view-layer';

const hazardLayerState = selector<ViewLayer[]>({
  key: 'hazardLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('hazards'))
      ? truthyKeys(get(hazardVisibilityState)).map((hazard) =>
          hazardViewLayer(hazard, get(dataParamsByGroupState(hazard)) as HazardParams),
        )
      : [],
});

const networkLayersState = selector<ViewLayer[]>({
  key: 'networkLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('assets'))
      ? get(networkSelectionState).map((network) => INFRASTRUCTURE_VIEW_LAYERS[network])
      : [],
});

export const viewLayersState = selector<ConfigTree<ViewLayer>>({
  key: 'viewLayersState',
  get: ({ get }) => {
    const showRegions = get(sectionVisibilityState('regions'));
    const regionLevel = get(regionLevelState);
    const background = get(backgroundState);
    const showLabels = get(showLabelsState);
    const isRetina = get(isRetinaState);

    return [
      // administrative region boundaries or population density
      showRegions &&
        (get(showPopulationState) ? populationViewLayer(regionLevel) : regionBoundariesViewLayer(regionLevel)),
      // hazard data layers
      get(hazardLayerState),

      get(sectionVisibilityState('buildings')) && buildingsViewLayer(),

      // network data layers
      get(networkLayersState),

      showLabels && [
        // basemap labels
        viewOnlyLayer('labels', () => labelsLayer(isRetina)),

        // administrative regions labels
        showRegions &&
          viewOnlyLayer(`boundaries_${regionLevel}-text`, () => regionLabelsDeckLayer(regionLevel, background)),
      ],
    ];
  },
});
