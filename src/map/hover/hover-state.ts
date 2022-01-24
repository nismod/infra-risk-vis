import { atom } from 'recoil';
import { RegionHover } from '../DataMap';

export const regionHoverState = atom<RegionHover>({
  key: 'regionHover',
  default: null,
});
