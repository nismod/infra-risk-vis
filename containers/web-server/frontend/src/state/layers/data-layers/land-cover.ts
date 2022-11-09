import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { landCoverViewLayer } from '@/config/land-cover/land-cover-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';

export const landCoverLayerState = selector<ViewLayer>({
  key: 'landCoverLayerState',
  get: ({ get }) => get(sidebarPathVisibilityState('exposure/land-cover')) && landCoverViewLayer(),
});
