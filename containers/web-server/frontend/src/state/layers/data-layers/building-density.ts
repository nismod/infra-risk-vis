import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { buildingDensityLayer } from '@/config/building-density/building-density-layer';
import { buildingDensityTypeState } from '@/state/data-selection/building-density';
import { sectionVisibilityState } from '@/state/sections';

export const buildingDensityLayerState = selector<ViewLayer>({
  key: 'buildingDensityLayerState',
  get: ({ get }) => get(sectionVisibilityState('buildings')) && buildingDensityLayer(get(buildingDensityTypeState)),
});
