import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { humanDevelopmentLayer } from '@/config/human-development/human-development-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { hdiRegionLevelState, hdiVariableState } from '@/state/data-selection/human-development';

export const humanDevelopmentLayerState = selector<ViewLayer>({
  key: 'humanDevelopmentLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('vulnerability/human/human-development')) &&
    humanDevelopmentLayer(get(hdiRegionLevelState), get(hdiVariableState)),
});
