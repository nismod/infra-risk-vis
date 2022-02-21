import { selector } from 'recoil';

import { flattenLayers } from 'lib/layer-tree';
import { ViewLayer } from 'lib/data-map/view-layers';

import { viewLayersState } from './view-layers';

export const viewLayersFlatState = selector<ViewLayer[]>({
  key: 'viewLayersFlatState',
  get: ({ get }) => flattenLayers(get(viewLayersState)),
});
