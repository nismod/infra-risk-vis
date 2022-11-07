import produce from 'immer';
import _ from 'lodash';
import { atom, selector } from 'recoil';

import { StateEffect } from '@/lib/recoil/state-effects/types';

import { HAZARD_DOMAINS_CONFIG } from '@/config/hazards/domains';
import { sidebarVisibilityToggleState } from '@/sidebar/SidebarContent';
import { paramsState } from '@/state/data-params';
import { viewState } from '@/state/view';

export const showDamagesState = selector({
  key: 'showDamagesState',
  get: ({ get }) => get(viewState) === 'risk', // get(networksStyleState) === 'damages',
});

export const damageSourceState = atom({
  key: 'damageSourceState',
  default: 'fluvial',
});

export const damageTypeState = atom({
  key: 'damageTypeState',
  default: 'direct',
});

const syncHazardsWithDamageSourceStateEffect = ({ get, set }, damageSource: string) => {
  _.forEach(HAZARD_DOMAINS_CONFIG, (groupConfig, group) => {
    set(sidebarVisibilityToggleState(`hazards/${group}`), group === damageSource);
  });
};

export const damageSourceStateEffect: StateEffect<string> = ({ get, set }, damageSource) => {
  syncHazardsWithDamageSourceStateEffect({ get, set }, damageSource);

  if (damageSource !== 'all') {
    const state = get(paramsState(damageSource));

    const damageSourceReturnPeriodDomain = state['rp'].options;
    const topReturnPeriod = damageSourceReturnPeriodDomain[damageSourceReturnPeriodDomain.length - 1];

    // CAUTION: this won't resolve the dependencies between data params if any depend on the return period
    set(
      paramsState(damageSource),
      produce(state, (draft) => {
        draft['rp'].value = topReturnPeriod;
      }),
    );
  }
};
