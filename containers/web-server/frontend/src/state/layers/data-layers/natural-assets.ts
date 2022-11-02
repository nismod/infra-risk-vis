import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { natureRasterViewLayer } from '@/config/natural-assets/nature-raster-layer';
import { naturalAssetsSelectionState } from '@/state/data-selection/natural-assets';
import { sectionVisibilityState } from '@/state/sections';

export const naturalAssetsLayersState = selector<ViewLayer[]>({
  key: 'naturalAssetsLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('natural-assets'))
      ? truthyKeys(get(naturalAssetsSelectionState)).map((type) => natureRasterViewLayer(type))
      : [],
});
