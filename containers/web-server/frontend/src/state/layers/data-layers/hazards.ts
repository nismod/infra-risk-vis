import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { hazardViewLayer } from '@/config/hazards/hazard-view-layer';
import { dataParamsByGroupState } from '@/state/data-params';
import { hazardVisibilityState } from '@/state/hazards/hazard-visibility';
import { sectionVisibilityState } from '@/state/sections';

export const hazardLayerState = selector<ViewLayer[]>({
  key: 'hazardLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('hazards'))
      ? truthyKeys(get(hazardVisibilityState)).map((hazard) =>
          hazardViewLayer(hazard, get(dataParamsByGroupState(hazard))),
        )
      : [],
});
