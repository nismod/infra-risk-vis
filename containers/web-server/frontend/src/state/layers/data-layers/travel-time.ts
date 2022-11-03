import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { travelTimeViewLayer } from '@/config/travel-time/travel-time-layer';
import { travelTimeTypeState } from '@/state/data-selection/travel-time';
import { sectionVisibilityState } from '@/state/sections';

export const travelTimeLayerState = selector<ViewLayer>({
  key: 'travelTimeLayerState',
  get: ({ get }) => get(sectionVisibilityState('travel-time')) && travelTimeViewLayer(get(travelTimeTypeState)),
});
