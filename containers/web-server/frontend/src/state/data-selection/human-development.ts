import { atom } from 'recoil';

import { HdiRegionLevel, HdiVariableType } from '@/config/human-development/metadata';

export const hdiVariableState = atom<HdiVariableType>({
  key: 'hdiVariableState',
  default: 'subnational_hdi',
});

export const hdiRegionLevelState = atom<HdiRegionLevel>({
  key: 'hdiRegionLevelState',
  default: 'countries',
});
