import { atom } from 'recoil';
import { BoundaryLevel } from '../../config/deck-layers/boundaries-layer';

import { BackgroundName } from '../../config/backgrounds';

export const backgroundState = atom<BackgroundName>({
  key: 'background',
  default: 'light',
});

export const showLabelsState = atom<boolean>({
  key: 'showLabels',
  default: true,
});

export const showBoundariesState = atom<boolean>({
  key: 'showBoundaries',
  default: false,
});

export const boundaryLevelState = atom<BoundaryLevel>({
  key: 'boundaryLevel',
  default: 'parish',
});
