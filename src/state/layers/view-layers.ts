import { regionLabelsDeckLayer } from 'config/regions/region-labels-deck-layer';
import { regionBoundariesViewLayer } from 'config/regions/boundaries-view-layer';
import { ViewLayer, viewOnlyLayer } from 'lib/data-map/view-layers';
import { backgroundState, showLabelsState } from 'map/layers/layers-state';
import { selector } from 'recoil';
import { truthyKeys } from 'lib/helpers';
import { labelsLayer } from 'config/deck-layers/labels-layer';
import { isRetinaState } from 'state/is-retina';
import { ConfigTree } from 'lib/nested-config/config-tree';

import { populationViewLayer } from 'config/regions/population-view-layer';
import { regionLevelState, showPopulationState } from 'state/regions';
import { sectionVisibilityState } from 'state/sections';
import { buildingsViewLayer } from 'config/buildings/buildings-view-layer';
import { buildingSelectionState } from 'state/buildings';
import { networkLayersState } from './networks';
import { hazardLayerState } from './hazards';

const buildingLayersState = selector<ViewLayer[]>({
  key: 'buildingLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('buildings'))
      ? truthyKeys(get(buildingSelectionState)).map((buildingType) => buildingsViewLayer(buildingType))
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

      get(buildingLayersState),

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
