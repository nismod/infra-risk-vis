import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { jrcPopulationViewLayer } from '@/config/population/population-view-layer';
import { sectionVisibilityState } from '@/state/sections';

export const populationLayerState = selector<ViewLayer>({
  key: 'populationLayerState',
  get: ({ get }) => get(sectionVisibilityState('population')) && jrcPopulationViewLayer(),
});
