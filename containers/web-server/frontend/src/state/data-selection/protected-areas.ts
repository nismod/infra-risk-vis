import { atom } from 'recoil';

import { ProtectedAreaType } from '@/config/protected-areas/metadata';

export const protectedAreaTypeSelectionState = atom<Record<ProtectedAreaType, boolean>>({
  key: 'protectedAreaTypeSelectionState',
  default: {
    land: true,
    marine: true,
  },
});
