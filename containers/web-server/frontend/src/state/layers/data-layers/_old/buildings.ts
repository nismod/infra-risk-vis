import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { buildingsViewLayer } from '@/config/_old/buildings/buildings-view-layer';
import { buildingSelectionState } from '@/state/data-selection/_old/buildings';
import { sectionVisibilityState } from '@/state/sections';

export const buildingLayersState = selector<ViewLayer[]>({
  key: 'buildingLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('buildings'))
      ? truthyKeys(get(buildingSelectionState)).map((buildingType) => buildingsViewLayer(buildingType))
      : [],
});
