import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { natureRasterViewLayer } from '@/config/natural-assets/nature-raster-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';

export const organicCarbonLayerState = selector<ViewLayer>({
  key: 'organicCarbonLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('exposure/organic-carbon')) ? natureRasterViewLayer('organic_carbon') : null,
});
