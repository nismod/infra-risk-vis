import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { healthsitesViewLayer } from '@/config/healthcare/healthsites-view-layer';
import { sectionVisibilityState } from '@/state/sections';

export const healthcareLayersState = selector<ViewLayer>({
  key: 'healthcareLayersState',
  get: ({ get }) => {
    return get(sectionVisibilityState('healthcare')) ? healthsitesViewLayer() : null;
  },
});
