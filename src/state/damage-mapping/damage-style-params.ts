import { damageSourceState } from './damage-map';
import { dataParamsByGroupState } from '../data-params';
import { selector } from 'recoil';
import { Accessor, withTriggers } from 'lib/deck/props/getters';

function getEadKey(hazard: string, rcp: string, epoch: number) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

/**
 * Returns sum of array elements, or null if all elements are null
 */
function sumOrNone(arr: number[]): number | null {
  let result: number = null;

  for (const x of arr) {
    if (x != null) {
      if (result == null) result = x;
      else result += x;
    }
  }
  return result;
}

function getEadAccessor(eadSource: string, rcp: string, epoch: number): Accessor<any> {
  let fn: any;
  if (eadSource === 'total-damages') {
    fn = (f) => {
      const values = ['fluvial', 'surface', 'coastal', 'cyclone'].map((ht) => f.properties[getEadKey(ht, rcp, epoch)]);
      return sumOrNone(values);
    };
  } else {
    fn = (f) => f.properties[getEadKey(eadSource, rcp, epoch)];
  }

  return withTriggers(fn, [eadSource, rcp, epoch]);
}

export const eadFieldSpecState = selector({
  key: 'eadAccessorState',
  get: ({ get }) => {
    const eadSource = get(damageSourceState);
    if (eadSource == null) return null;
    const eadParams = get(dataParamsByGroupState(eadSource));

    return {
      field: 'damages',
      fieldParams: {
        damage_type: 'direct',
        hazard: eadSource,
        rcp: eadParams.rcp,
        epoch: eadParams.epoch,
        protection_standard: 0,
      },
    };
    // return getEadAccessor(eadSource, eadParams.rcp, eadParams.epoch);
  },
});

export const damageMapStyleParamsState = selector({
  key: 'damageMapStyleParamsState',
  get: ({ get }) => {
    const eadFieldSpec = get(eadFieldSpecState);
    if (eadFieldSpec == null) return {};

    return {
      colorMap: {
        colorScheme: 'damages',
        colorField: eadFieldSpec,
      },
    };
  },
});
