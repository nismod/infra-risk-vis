import { selector } from 'recoil';

import { StyleParams, ViewLayer } from '@/lib/data-map/view-layers';

import { infrastructureViewLayer } from '@/config/networks/infrastructure-view-layer';
import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { damageMapStyleParamsState } from '@/state/data-selection/damage-mapping/damage-style-params';
import { adaptationStyleParamsState } from '@/state/data-selection/networks/adaptations';
import { networkSelectionState } from '@/state/data-selection/networks/network-selection';
import { networksStyleState } from '@/state/data-selection/networks/networks-style';

export const networkLayersState = selector<ViewLayer[]>({
  key: 'networkLayersState',
  get: ({ get }) =>
    get(sidebarPathVisibilityState('exposure/infrastructure'))
      ? get(networkSelectionState).map((network) => infrastructureViewLayer(network, get(networkStyleParamsState)))
      : [],
});

export const networkStyleParamsState = selector<StyleParams>({
  key: 'networkStyleParamsState',
  get: ({ get }) => {
    switch (get(networksStyleState)) {
      case 'damages':
        return get(damageMapStyleParamsState);
      case 'adaptation':
        return get(adaptationStyleParamsState);
      default:
        return {};
    }
  },
});
