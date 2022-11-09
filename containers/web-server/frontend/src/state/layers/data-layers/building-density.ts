import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { buildingDensityLayer } from '@/config/building-density/building-density-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { buildingDensityTypeState } from '@/state/data-selection/building-density';

export const buildingDensityLayerState = selector<ViewLayer>({
  key: 'buildingDensityLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('exposure/buildings')) && buildingDensityLayer(get(buildingDensityTypeState)),
});
