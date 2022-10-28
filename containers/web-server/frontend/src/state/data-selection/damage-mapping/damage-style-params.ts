import { selector } from 'recoil';

import { FieldSpec, StyleParams } from '@/lib/data-map/view-layers';

import { DAMAGE_COLORMAP } from '@/config/damage-mapping/colors';
import { dataParamsByGroupState } from '@/state/data-params';

import { damageSourceState, damageTypeState } from './damage-map';

export const damagesFieldState = selector<FieldSpec>({
  key: 'eadAccessorState',
  get: ({ get }) => {
    const damageSource = get(damageSourceState);
    if (damageSource == null) return null;
    const damageType = get(damageTypeState);
    const damageParams = get(dataParamsByGroupState(damageSource));

    return {
      fieldGroup: 'damages_expected',
      fieldDimensions: {
        hazard: damageSource,
        rcp: damageParams.rcp,
        epoch: damageParams.epoch,
        protection_standard: 0,
      },
      field: damageType === 'direct' ? 'ead_mean' : 'eael_mean',
    };
  },
});

export const damageMapStyleParamsState = selector<StyleParams>({
  key: 'damageMapStyleParamsState',
  get: ({ get }) => {
    const eadFieldSpec = get(damagesFieldState);
    if (eadFieldSpec == null) return {};

    return {
      colorMap: {
        colorSpec: DAMAGE_COLORMAP,
        fieldSpec: eadFieldSpec,
      },
    };
  },
});
