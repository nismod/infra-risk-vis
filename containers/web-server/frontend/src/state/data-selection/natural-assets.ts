import { atom } from 'recoil';

export type NaturalAssetType = 'organic_carbon';

export const naturalAssetsSelectionState = atom<Record<NaturalAssetType, boolean>>({
  key: 'naturalAssetsSelectionState',
  default: {
    organic_carbon: false,
  },
});
