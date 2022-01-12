import { atom } from 'recoil';

export const placeSearchActiveState = atom<boolean>({
  key: 'placeSearchActive',
  default: false,
});

export const placeSearchQueryState = atom<string>({
  key: 'placeSearchQuery',
  default: '',
});

export interface PlaceSearchResult {
  label: string;
  latitude: number;
  longitude: number;
}

export const placeSearchSelectedResultState = atom<PlaceSearchResult>({
  key: 'placeSearchSelectedResult',
  default: null,
});
