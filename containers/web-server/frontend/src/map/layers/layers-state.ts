import { atom } from 'recoil';

import { BackgroundName } from '@/config/backgrounds';

export const backgroundState = atom<BackgroundName>({
  key: 'background',
  default: 'light',
});

export const showLabelsState = atom<boolean>({
  key: 'showLabels',
  default: true,
});
