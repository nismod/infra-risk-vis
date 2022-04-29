import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { AdaptationOptionParams } from 'config/domains/adaptation';
import { INFRASTRUCTURE_VIEW_LAYERS } from 'config/networks/view-layers';
import { ViewLayer, StyleParams, ColorSpec } from 'lib/data-map/view-layers';
import { StateEffect } from 'lib/recoil/state-effects/types';
import { atom, selector } from 'recoil';
import { damageMapStyleParamsState } from 'state/damage-mapping/damage-style-params';
import { dataParamsByGroupState } from 'state/data-params';
import { networkSelectionState, networkTreeCheckboxState, networkTreeConfig } from 'state/networks/network-selection';
import { networksStyleState } from 'state/networks/networks-style';
import { sectionVisibilityState } from 'state/sections';

import adaptationSectorLayers from 'config/domains/adaptation-sector-layers.json';
import _ from 'lodash';
import { recalculateCheckboxStates } from 'lib/controls/checkbox-tree/CheckboxTree';

export const networkLayersState = selector<ViewLayer[]>({
  key: 'networkLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('assets'))
      ? get(networkSelectionState).map((network) => INFRASTRUCTURE_VIEW_LAYERS[network])
      : [],
});

export const showAdaptationsState = selector<boolean>({
  key: 'showAdaptationsState',
  get: ({ get }) => get(networksStyleState) === 'adaptation',
});

export const adaptationFieldState = atom<'avoided_ead_mean' | 'avoided_eael_mean' | 'adaptation_cost'>({
  key: 'adaptationFieldState',
  default: 'avoided_ead_mean',
});

export const adaptationDataParamsStateEffect: StateEffect<AdaptationOptionParams> = (
  { get, set },
  adaptationParams,
) => {
  const { sector, subsector, asset_type } = adaptationParams;

  const layers = _.uniq(
    adaptationSectorLayers
      .filter((x) => x.sector === sector && x.subsector === subsector && x.asset_type === asset_type)
      .map((x) => x.layer_name),
  );

  const currentSelection = get(networkTreeCheckboxState);
  const updatedTreeState = {
    checked: {
      ..._.mapValues(currentSelection.checked, () => false),
      ..._.fromPairs(layers.map((layer) => [layer, true])),
    },
    indeterminate: {},
  };
  const resolvedTreeState = recalculateCheckboxStates(updatedTreeState, networkTreeConfig);

  set(networkTreeCheckboxState, resolvedTreeState);

  // currently not auto-updating the expanded state of the tree since that can make the adaptations UI section move out of view
  // set(networkTreeExpandedState, truthyKeys(resolvedTreeState.indeterminate));
};

export const adaptationStyleParamsState = selector<StyleParams>({
  key: 'adaptationStyleParamsState',
  get: ({ get }) => {
    const field = get(adaptationFieldState);
    const { hazard, rcp, adaptation_name, adaptation_protection_level } = get(dataParamsByGroupState('adaptation'));

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
            hazard,
            rcp,
            adaptation_name,
            adaptation_protection_level,
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
