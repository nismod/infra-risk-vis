import _ from 'lodash';
import { TransactionInterface_UNSTABLE, atom, selector } from 'recoil';

import { recalculateCheckboxStates } from '@/lib/controls/checkbox-tree/CheckboxTree';
import { d3Scale, d3ScaleChromatic, discardSides, invertColorScale } from '@/lib/data-map/color-maps';
import { ColorSpec, FieldSpec, StyleParams } from '@/lib/data-map/view-layers';
import { valueType } from '@/lib/helpers';
import { StateEffect } from '@/lib/recoil/state-effects/types';

import { LayerSpec } from '@/asset-list/use-sorted-features';
import { AdaptationOptionParams } from '@/config/domains/adaptation';
import adaptationSectorLayers from '@/config/domains/adaptation-sector-layers.json';
import { dataParamsByGroupState } from '@/state/data-params';

import { networkTreeCheckboxState, networkTreeConfig } from './network-selection';
import { networksStyleState } from './networks-style';

export const showAdaptationsState = selector<boolean>({
  key: 'showAdaptationsState',
  get: ({ get }) => get(networksStyleState) === 'adaptation',
});

export const adaptationFieldState = atom<
  'avoided_ead_mean' | 'avoided_eael_mean' | 'adaptation_cost' | 'cost_benefit_ratio'
>({
  key: 'adaptationFieldState',
  default: 'avoided_ead_mean',
});

export const adaptationCostBenefitRatioEaelDaysState = atom<number>({
  key: 'adaptationCostBenefitRatioEaelDaysState',
  default: 15,
});

export function syncInfrastructureSelectionStateEffect({ get, set }: TransactionInterface_UNSTABLE, layers: string[]) {
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
}

export const adaptationDataParamsStateEffect: StateEffect<AdaptationOptionParams> = (iface, adaptationParams) => {
  const { sector, subsector, asset_type } = adaptationParams;

  const layers = _.uniq(
    adaptationSectorLayers
      .filter((x) => x.sector === sector && x.subsector === subsector && x.asset_type === asset_type)
      .map((x) => x.layer_name),
  );

  syncInfrastructureSelectionStateEffect(iface, layers);

  // currently not auto-updating the expanded state of the tree since that can make the adaptations UI section move out of view
  // set(networkTreeExpandedState, truthyKeys(resolvedTreeState.indeterminate));
};

export const adaptationLayerSpecState = selector<LayerSpec>({
  key: 'adaptationLayerSpecState',
  get: ({ get }) => {
    const { sector, subsector, asset_type } = get(dataParamsByGroupState('adaptation'));

    return {
      sector,
      subsector,
      assetType: asset_type,
    };
  },
});

export const adaptationFieldSpecState = selector<FieldSpec>({
  key: 'adaptationFieldSpecState',
  get: ({ get }) => {
    const field = get(adaptationFieldState);
    const { hazard, rcp, adaptation_name, adaptation_protection_level } = get(dataParamsByGroupState('adaptation'));

    let fieldParams: any = {};
    if (field === 'cost_benefit_ratio') {
      fieldParams = {
        eael_days: get(adaptationCostBenefitRatioEaelDaysState),
      };
    }

    return {
      fieldGroup: 'adaptation',
      fieldDimensions: {
        hazard,
        rcp,
        adaptation_name,
        adaptation_protection_level,
      },
      field,
      fieldParams,
    };
  },
});

const ADAPTATION_COLOR_MAPS = valueType<ColorSpec>()({
  cost: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateGreens, 0.2, 0.2),
    range: [0, 1e9],
    empty: '#ccc',
  },
  avoided: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 1e7],
    empty: '#ccc',
  },
  costBenefitRatio: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateViridis),
    range: [1, 10],
    empty: '#ccc',
  },
});
export const adaptationColorSpecState = selector<ColorSpec>({
  key: 'adaptationColorSpecState',
  get: ({ get }) => {
    const field = get(adaptationFieldState);

    let colorSpec: ColorSpec;
    if (field === 'adaptation_cost') {
      colorSpec = ADAPTATION_COLOR_MAPS.cost;
    } else if (field === 'avoided_ead_mean' || field === 'avoided_eael_mean') {
      colorSpec = ADAPTATION_COLOR_MAPS.avoided;
    } else if (field === 'cost_benefit_ratio') {
      colorSpec = ADAPTATION_COLOR_MAPS.costBenefitRatio;
    }

    return colorSpec;
  },
});

export const adaptationStyleParamsState = selector<StyleParams>({
  key: 'adaptationStyleParamsState',
  get: ({ get }) => {
    const fieldSpec = get(adaptationFieldSpecState);
    const colorSpec = get(adaptationColorSpecState);

    return {
      colorMap: {
        fieldSpec,
        colorSpec,
      },
    };
  },
});
