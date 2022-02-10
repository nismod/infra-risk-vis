import { atom } from 'recoil';

import { PlaceSearchResult } from 'lib/map/place-search/use-place-search';

export const placeSearchActiveState = atom<boolean>({
  key: 'placeSearchActive',
  default: false,
});

export const placeSearchQueryState = atom<string>({
  key: 'placeSearchQuery',
  default: '',
});

export const placeSearchSelectedResultState = atom<PlaceSearchResult>({
  key: 'placeSearchSelectedResult',
  default: null,
});
