import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { travelTimeViewLayer } from '@/config/travel-time/travel-time-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { travelTimeTypeState } from '@/state/data-selection/travel-time';

export const travelTimeLayerState = selector<ViewLayer>({
  key: 'travelTimeLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('vulnerability/human/travel-time')) && travelTimeViewLayer(get(travelTimeTypeState)),
});
