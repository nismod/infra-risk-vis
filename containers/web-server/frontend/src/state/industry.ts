import { atom } from 'recoil';

export const INDUSTRY_TYPES = ['cement', 'steel'] as const;

export type IndustryType = typeof INDUSTRY_TYPES[number];

export type IndustrySelection = Record<IndustryType, boolean>;

export const industrySelectionState = atom<IndustrySelection>({
  key: 'industrySelectionState',
  default: {
    cement: true,
    steel: true,
  },
});
