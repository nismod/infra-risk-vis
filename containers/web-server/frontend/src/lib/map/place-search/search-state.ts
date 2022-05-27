import { atom } from 'recoil';

export const placeSearchActiveState = atom<boolean>({
  key: 'placeSearchActive',
  default: false,
});

export const placeSearchQueryState = atom<string>({
  key: 'placeSearchQuery',
  default: '',
});
