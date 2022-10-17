import { selector } from 'recoil';

import { INTERACTION_GROUPS } from '@/config/interaction-groups';
import { showPopulationState } from '@/state/regions';

export const interactionGroupsState = selector({
  key: 'interactionGroupsState',
  get: ({ get }) => {
    const regionDataShown = get(showPopulationState);

    return [
      INTERACTION_GROUPS.assets,
      INTERACTION_GROUPS.hazards,
      {
        ...INTERACTION_GROUPS.regions,
        usesAutoHighlight: regionDataShown,
      },
      INTERACTION_GROUPS.solutions,
      INTERACTION_GROUPS.drought,
    ];
  },
});
