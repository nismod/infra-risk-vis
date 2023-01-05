import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { regionalExposureLayer } from '@/config/regional-risk/regional-risk-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { regionalExposureVariableState } from '@/state/data-selection/regional-risk';

export const regionalExposureLayerState = selector<ViewLayer>({
  key: 'regionalExposureLayerState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('risk/regional')) &&
    regionalExposureLayer(get(regionalExposureVariableState)),
});
