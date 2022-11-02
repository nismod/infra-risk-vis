import { atom } from 'recoil';

import { TraveltimeType } from '@/config/travel-time/travel-time-layer';

export const travelTimeTypeState = atom<TraveltimeType>({
  key: 'travelTimeTypeState',
  default: 'motorized',
});
