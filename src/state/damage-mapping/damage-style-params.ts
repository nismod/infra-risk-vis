import { damageSourceState, damageTypeState } from './damage-map';
import { dataParamsByGroupState } from '../data-params';
import { selector } from 'recoil';
import { FieldSpec } from 'lib/data-map/view-layers';

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

export const damageMapStyleParamsState = selector({
  key: 'damageMapStyleParamsState',
  get: ({ get }) => {
    const eadFieldSpec = get(damagesFieldState);
    if (eadFieldSpec == null) return {};

    return {
      colorMap: {
        colorScheme: 'damages',
        colorField: eadFieldSpec,
      },
    };
  },
});
