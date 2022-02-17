import _ from 'lodash';
import { selector } from 'recoil';

import { selectedDamageSourceState, showDirectDamagesState } from './damage-mapping/damage-map';
import { dataParamsByGroupState } from './data-params';

function getEadKey(hazard: string, rcp: string, epoch: number) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

function getEadAccessor(eadSource: string, rcp: string, epoch: number) {
  if (eadSource === 'total-damages') {
    return (f) =>
      _.sum(['fluvial', 'surface', 'coastal', 'cyclone'].map((ht) => f.properties[getEadKey(ht, rcp, epoch)] ?? 0));
  } else {
    return getEadKey(eadSource, rcp, epoch);
  }
}

function getDamageMapStyleParams(eadSource, eadSourceParams) {
  const { rcp, epoch } = eadSourceParams;

  return {
    colorMap: {
      colorScheme: 'damages',
      colorField: getEadAccessor(eadSource, rcp, epoch),
    },
  };
}

const damageMapStyleParamsState = selector({
  key: 'damageMapStyleParamsState',
  get: ({ get }) => {
    const eadSource = get(selectedDamageSourceState);
    if (eadSource == null) return {};
    const eadSourceParams = get(dataParamsByGroupState(eadSource));

    return getDamageMapStyleParams(eadSource, eadSourceParams);
  },
});

export const styleParamsState = selector({
  key: 'styleParamsState',
  get: ({ get }) => {
    if (!get(showDirectDamagesState)) return {};

    return get(damageMapStyleParamsState);
  },
});
