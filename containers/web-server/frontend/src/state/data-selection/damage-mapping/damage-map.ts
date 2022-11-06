import _ from 'lodash';
import { atom, selector } from 'recoil';

import { HAZARD_DOMAINS_CONFIG } from '@/config/hazards/domains';
import { paramOptionsState, paramValueState } from '@/state/data-params';
import { hazardToggleState } from '@/state/data-selection/hazards/hazard-selection';
import { networksStyleState } from '@/state/data-selection/networks/networks-style';

export const showDamagesState = selector({
  key: 'showDamagesState',
  get: ({ get }) => get(networksStyleState) === 'damages',
});

export const damageSourceState = atom({
  key: 'damageSourceState',
  default: 'all',
});

export const damageTypeState = atom({
  key: 'damageTypeState',
  default: 'direct',
});

export const damageSourceStateEffect = ({ get, set }, damageSource) => {
  syncHazardsWithDamageSourceStateEffect({ get, set }, damageSource);

  if (damageSource !== 'all') {
    const damageSourceReturnPeriodDomain = get(paramOptionsState({ group: damageSource, param: 'rp' }));
    const topReturnPeriod = damageSourceReturnPeriodDomain[damageSourceReturnPeriodDomain.length - 1];

    // CAUTION: this won't resolve the dependencies between data params if any depend on the return period
    set(paramValueState({ group: damageSource, param: 'rp' }), topReturnPeriod);
  }
};

function syncHazardsWithDamageSourceStateEffect({ get, set }, damageSource) {
  _.forEach(HAZARD_DOMAINS_CONFIG, (groupConfig, group) => {
    set(hazardToggleState(group), group === damageSource);
  });
}
