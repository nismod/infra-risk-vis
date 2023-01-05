import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { jrcPopulationViewLayer } from '@/config/population/population-view-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';

export const populationLayerState = selector<ViewLayer>({
  key: 'populationLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('exposure/population')) && jrcPopulationViewLayer(),
});
