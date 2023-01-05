import produce from 'immer';
import _ from 'lodash';
import { atom, selector } from 'recoil';

import { CurrentStateEffect, StateEffect } from '@/lib/recoil/state-effects/types';

import { HAZARD_DOMAINS_CONFIG } from '@/config/hazards/domains';
import { sidebarVisibilityToggleState } from '@/sidebar/SidebarContent';
import { paramsState } from '@/state/data-params';
import { viewState } from '@/state/view';

export const showInfrastructureDamagesState = selector({
  key: 'showInfrastructureDamagesState',
  get: ({ get }) =>
    get(viewState) === 'risk' && get(sidebarVisibilityToggleState('risk/infrastructure')),
});

export const damageSourceState = atom({
  key: 'damageSourceState',
  default: 'fluvial',
});

export const damageTypeState = atom({
  key: 'damageTypeState',
  default: 'direct',
});

export const syncHazardsWithDamageSourceStateEffect: CurrentStateEffect<string> = (
  { get, set },
  damageSource,
) => {
  _.forEach(HAZARD_DOMAINS_CONFIG, (groupConfig, group) => {
    set(sidebarVisibilityToggleState(`hazards/${group}`), group === damageSource);
  });
};

export const damageSourceStateEffect: StateEffect<string> = (iface, damageSource) => {
  syncHazardsWithDamageSourceStateEffect(iface, damageSource);
  const { get, set } = iface;

  if (damageSource !== 'all') {
    /**
     * If damage source is not 'all', set return period for hazard to highest available.
     */
    const state = get(paramsState(damageSource));

    const damageSourceReturnPeriodDomain = state['rp'].options;
    const topReturnPeriod =
      damageSourceReturnPeriodDomain[damageSourceReturnPeriodDomain.length - 1];

    // CAUTION: this won't resolve the dependencies between data params if any depend on the return period
    set(
      paramsState(damageSource),
      produce(state, (draft) => {
        draft['rp'].value = topReturnPeriod;
      }),
    );
  }
};
