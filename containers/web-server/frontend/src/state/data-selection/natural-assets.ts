import { atom } from 'recoil';

import { NaturalAssetType } from '@/config/natural-assets/metadata';

export const naturalAssetsSelectionState = atom<Record<NaturalAssetType, boolean>>({
  key: 'naturalAssetsSelectionState',
  default: {
    biodiversity_intactness: true,
    forest_landscape_integrity: false,
    organic_carbon: false,
  },
});
