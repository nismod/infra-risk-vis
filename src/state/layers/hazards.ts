import { HazardParams } from 'config/hazards/domains';
import { hazardViewLayer } from 'config/hazards/hazard-view-layer';
import { ViewLayer } from 'lib/data-map/view-layers';
import { truthyKeys } from 'lib/helpers';
import { selector } from 'recoil';
import { dataParamsByGroupState } from 'state/data-params';
import { hazardVisibilityState } from 'state/hazards/hazard-visibility';
import { sectionVisibilityState } from 'state/sections';

export const hazardLayerState = selector<ViewLayer[]>({
  key: 'hazardLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('hazards'))
      ? truthyKeys(get(hazardVisibilityState)).map((hazard) =>
          hazardViewLayer(hazard, get(dataParamsByGroupState(hazard)) as HazardParams),
        )
      : [],
});
