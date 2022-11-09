import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { industryViewLayer } from '@/config/industry/industry-view-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { industrySelectionState } from '@/state/data-selection/industry';

export const industryLayersState = selector<ViewLayer[]>({
  key: 'industryLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('exposure/industry'))
      ? truthyKeys(get(industrySelectionState)).map((industryType) => industryViewLayer(industryType))
      : [],
});
