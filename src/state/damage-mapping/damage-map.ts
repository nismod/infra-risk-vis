import { HAZARD_DOMAINS } from 'config/hazards/domains';
import _ from 'lodash';
import { atom, selector, useRecoilTransaction_UNSTABLE } from 'recoil';
import { dataParamOptionsState, dataParamState } from 'state/data-params';
import { hazardSelectionState } from 'state/hazards/hazard-selection';
import { networksStyleState } from 'state/networks/networks-style';

export const showDirectDamagesState = selector({
  key: 'showDamagesState',
  get: ({ get }) => get(networksStyleState) === 'direct-damages',
});

const selectedDamageSourceImpl = atom({
  key: 'selectedDamageSourceImpl',
  default: 'total-damages',
});

export const selectedDamageSourceState = selector<string>({
  key: 'selectedDamageSourceState',
  get: ({ get }) => get(selectedDamageSourceImpl),
});

function syncHazardsWithDamageSourceTask({ get, set }) {
  const damageSource = get(selectedDamageSourceImpl);

  _.forEach(HAZARD_DOMAINS, (groupConfig, group) => {
    set(hazardSelectionState(group), group === damageSource);
  });
}

export function useInitDirectDamages() {
  return useRecoilTransaction_UNSTABLE(({ get, set }) => () => {
    syncHazardsWithDamageSourceTask({ get, set });
  });
}

export function useUpdateDamageSource() {
  return useRecoilTransaction_UNSTABLE(({ get, set }) => (damageSource: string) => {
    // const currentDamageSource = get(selectedDamageSourceImpl);
    set(selectedDamageSourceImpl, damageSource);

    syncHazardsWithDamageSourceTask({ get, set });

    // if (currentDamageSource !== 'total-damages') {
    //   set(hazardSelectionState(currentDamageSource), false);
    // }
    if (damageSource !== 'total-damages') {
      // set(hazardSelectionState(damageSource), true);
      const damageSourceReturnPeriodDomain = get(dataParamOptionsState({ group: damageSource, param: 'returnPeriod' }));
      const topReturnPeriod = damageSourceReturnPeriodDomain[damageSourceReturnPeriodDomain.length - 1];

      // CAUTION: this won't resolve the dependencies between data params if any depend on the return period
      set(dataParamState({ group: damageSource, param: 'returnPeriod' }), topReturnPeriod);
    }
  });
}
