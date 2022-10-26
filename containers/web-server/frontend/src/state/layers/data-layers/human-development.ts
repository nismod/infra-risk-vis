import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { humanDevelopmentLayer } from '@/config/human-development/human-development-layer';
import { hdiRegionLevelState, hdiVariableState } from '@/state/data-selection/human-development';
import { sectionVisibilityState } from '@/state/sections';

export const humanDevelopmentLayerState = selector<ViewLayer>({
  key: 'humanDevelopmentLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('human-vulnerability')) &&
    humanDevelopmentLayer(get(hdiRegionLevelState), get(hdiVariableState)),
});
