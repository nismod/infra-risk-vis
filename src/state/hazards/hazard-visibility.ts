import { selector } from 'recoil';

import { HAZARDS_MAP_ORDER } from 'config/hazards/metadata';
import { damageSourceState, showDirectDamagesState } from 'state/damage-mapping/damage-map';
import { getHazardSelectionAggregate } from './hazard-selection';

export const hazardVisibilityState = selector({
  key: 'hazardVisibilityState',
  get: ({ get }) => {
    if (get(showDirectDamagesState)) {
      const selectedDamageSource = get(damageSourceState);
      if (selectedDamageSource === 'total-damages') {
        return {};
      } else {
        return {
          [selectedDamageSource]: true,
        };
      }
    } else {
      return getHazardSelectionAggregate({ get }, HAZARDS_MAP_ORDER);
    }
  },
});
