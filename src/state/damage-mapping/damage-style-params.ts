import { damageSourceState } from './damage-map';
import { dataParamsByGroupState } from '../data-params';
import { selector } from 'recoil';

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

function getEadAccessor(eadSource: string, rcp: string, epoch: number) {
  if (eadSource === 'total-damages') {
    return (f) => {
      const values = ['fluvial', 'surface', 'coastal', 'cyclone'].map((ht) => f.properties[getEadKey(ht, rcp, epoch)]);
      return sumOrNone(values);
    };
  } else {
    return (f) => f.properties[getEadKey(eadSource, rcp, epoch)];
  }
}

export const eadAccessorState = selector({
  key: 'eadAccessorState',
  get: ({ get }) => {
    const eadSource = get(damageSourceState);
    if (eadSource == null) return null;
    const eadParams = get(dataParamsByGroupState(eadSource));

    return getEadAccessor(eadSource, eadParams.rcp, eadParams.epoch);
  },
});

export const damageMapStyleParamsState = selector({
  key: 'damageMapStyleParamsState',
  get: ({ get }) => {
    const eadAccessor = get(eadAccessorState);
    if (eadAccessor == null) return {};

    return {
      colorMap: {
        colorScheme: 'damages',
        colorField: eadAccessor,
      },
    };
  },
});
