import { atom } from 'recoil';

import { BackgroundName } from '../../config/backgrounds';

export const backgroundState = atom<BackgroundName>({
  key: 'background',
  default: 'satellite',
});

export const showLabelsState = atom<boolean>({
  key: 'showLabels',
  default: true,
});
