import { selector } from 'recoil';

import { INTERACTION_GROUPS } from '@/config/interaction-groups';

export const interactionGroupsState = selector({
  key: 'interactionGroupsState',
  get: ({ get }) => {
    return [
      INTERACTION_GROUPS.assets,
      INTERACTION_GROUPS.wdpa,
      INTERACTION_GROUPS.hdi,
      INTERACTION_GROUPS.hazards,
      INTERACTION_GROUPS.population,
    ];
  },
});
