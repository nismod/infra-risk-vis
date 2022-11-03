import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { industryViewLayer } from '@/config/industry/industry-view-layer';
import { industrySelectionState } from '@/state/data-selection/industry';
import { sectionVisibilityState } from '@/state/sections';

export const industryLayersState = selector<ViewLayer[]>({
  key: 'industryLayerState',
  get: ({ get }) =>
    get(sectionVisibilityState('industry'))
      ? truthyKeys(get(industrySelectionState)).map((industryType) => industryViewLayer(industryType))
      : [],
});
