import { atom } from 'recoil';

import { BuildingDensityType } from '@/config/building-density/metadata';

export const buildingDensityTypeState = atom<BuildingDensityType>({
  key: 'buildingDensityTypeState',
  default: 'all',
});
