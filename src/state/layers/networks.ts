import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { INFRASTRUCTURE_VIEW_LAYERS } from 'config/networks/view-layers';
import { ViewLayer, StyleParams, ColorSpec } from 'lib/data-map/view-layers';
import { atom, selector } from 'recoil';
import { damageMapStyleParamsState } from 'state/damage-mapping/damage-style-params';
import { networkSelectionState } from 'state/networks/network-selection';
import { networksStyleState } from 'state/networks/networks-style';
import { sectionVisibilityState } from 'state/sections';

export const networkLayersState = selector<ViewLayer[]>({
  key: 'networkLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('assets'))
      ? get(networkSelectionState).map((network) => INFRASTRUCTURE_VIEW_LAYERS[network])
      : [],
});

export const adaptationFieldState = atom<'avoided_ead_mean' | 'avoided_eael_mean' | 'adaptation_cost'>({
  key: 'adaptationFieldState',
  default: 'avoided_ead_mean',
});

export const adaptationStyleParamsState = selector<StyleParams>({
  key: 'adaptationStyleParamsState',
  get: ({ get }) => {
    const field = get(adaptationFieldState);

    let colorSpec: ColorSpec;
    if (field === 'adaptation_cost') {
      colorSpec = VECTOR_COLOR_MAPS.adaptationCost;
    } else if (field === 'avoided_ead_mean' || field === 'avoided_eael_mean') {
      colorSpec = VECTOR_COLOR_MAPS.adaptationAvoided;
    }

    return {
      colorMap: {
        fieldSpec: {
          fieldGroup: 'adaptation',
          fieldDimensions: {
            hazard: 'flooding',
            rcp: '8.5',
            adaptation_name: 'Elevate the roads',
            adaptation_protection_level: 1,
          },
          field,
        },
        colorSpec,
      },
    };
  },
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
