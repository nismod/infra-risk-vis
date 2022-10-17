import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { flattenConfig } from '@/lib/nested-config/flatten-config';

import { viewLayersState } from './view-layers';

export const viewLayersFlatState = selector<ViewLayer[]>({
  key: 'viewLayersFlatState',
  get: ({ get }) => flattenConfig(get(viewLayersState)),
});
