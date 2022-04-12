import { damageSourceState, damageTypeState } from './damage-map';
import { dataParamsByGroupState } from '../data-params';
import { selector } from 'recoil';

export const damagesFieldState = selector({
  key: 'eadAccessorState',
  get: ({ get }) => {
    const damageSource = get(damageSourceState);
    if (damageSource == null) return null;
    const damageType = get(damageTypeState);
    const damageParams = get(dataParamsByGroupState(damageSource));

    return {
      field: 'damages',
      fieldParams: {
        damage_type: damageType,
        hazard: damageSource,
        rcp: damageParams.rcp,
        epoch: damageParams.epoch,
        protection_standard: 0,
      },
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
