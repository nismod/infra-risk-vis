import { selector } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { truthyKeys } from '@/lib/helpers';

import { ProtectedAreaType } from '@/config/protected-areas/metadata';
import { protectedAreaViewLayer } from '@/config/protected-areas/protected-area-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { protectedAreaTypeSelectionState } from '@/state/data-selection/protected-areas';

export const protectedAreasKeysState = selector<ProtectedAreaType[]>({
  key: 'protectedAreasKeysState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('vulnerability/nature/protected-areas'))
      ? truthyKeys(get(protectedAreaTypeSelectionState))
      : [],
});

export const protectedAreasPointLayerState = selector<ViewLayer[]>({
  key: 'protectedAreasPointLayerState',
  get: ({ get }) =>
    get(protectedAreasKeysState).map((type) => protectedAreaViewLayer('points', type)),
});

export const protectedAreasPolygonLayerState = selector<ViewLayer[]>({
  key: 'protectedAreasPolygonLayerState',
  get: ({ get }) =>
    get(protectedAreasKeysState).map((type) => protectedAreaViewLayer('polygons', type)),
});
