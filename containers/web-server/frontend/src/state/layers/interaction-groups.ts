import { selector } from 'recoil';

import { INTERACTION_GROUPS } from '@/config/interaction-groups';

export const interactionGroupsState = selector({
  key: 'interactionGroupsState',
  get: ({ get }) => {
    // the first group will be treated as primary and the picking radius from that group will be used globally
    return [
      INTERACTION_GROUPS.assets,
      INTERACTION_GROUPS.wdpa,
      INTERACTION_GROUPS.hdi,
      INTERACTION_GROUPS.hazards,
      INTERACTION_GROUPS.raster_assets,
    ];
  },
});
