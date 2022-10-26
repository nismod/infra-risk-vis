import { selector } from 'recoil';

import { INTERACTION_GROUPS } from '@/config/interaction-groups';

export const interactionGroupsState = selector({
  key: 'interactionGroupsState',
  get: ({ get }) => {
    return [
      INTERACTION_GROUPS.regions,
      INTERACTION_GROUPS.assets,
      INTERACTION_GROUPS.hazards,
      INTERACTION_GROUPS.population,
    ];
  },
});
