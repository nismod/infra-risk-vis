import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { healthsitesViewLayer } from '@/config/healthcare/healthsites-view-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';

export const healthcareLayersState = selector<ViewLayer>({
  key: 'healthcareLayersState',
  get: ({ get }) => {
    return get(sidebarPathVisibilityState('exposure/healthsites')) && healthsitesViewLayer();
  },
});
