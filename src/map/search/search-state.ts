import { atom } from 'recoil';

export const placeSearchActiveState = atom<boolean>({
  key: 'placeSearchActive',
  default: false,
});

export const placeSearchQueryState = atom<string>({
  key: 'placeSearchQuery',
  default: '',
});

export interface BoundingBox {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

export interface PlaceSearchResult {
  label: string;
  latitude: number;
  longitude: number;
  boundingBox: BoundingBox;
}

export const placeSearchSelectedResultState = atom<PlaceSearchResult>({
  key: 'placeSearchSelectedResult',
  default: null,
});
