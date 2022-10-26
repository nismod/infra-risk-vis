import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { protectedAreaViewLayer } from '@/config/protected-areas/protected-area-layer';
import { protectedAreaTypeSelectionState } from '@/state/data-selection/protected-areas';
import { sectionVisibilityState } from '@/state/sections';

export const protectedAreasLayerState = selector<ViewLayer[]>({
  key: 'protectedAreasLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('nature-vulnerability'))
      ? truthyKeys(get(protectedAreaTypeSelectionState)).map((type) => protectedAreaViewLayer(type))
      : [],
});
