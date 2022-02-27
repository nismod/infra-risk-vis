import { RegionLevel } from 'config/regions/metadata';
import { atom, selector } from 'recoil';

export const showRegionsState = atom<boolean>({
  key: 'showRegionsState',
  default: true,
});

export const regionLevelState = atom<RegionLevel>({
  key: 'regionLevelState',
  default: 'parish',
});

export type RegionDataType = 'boundaries' | 'population';

export const regionDataState = atom<RegionDataType>({
  key: 'regionDataState',
  default: 'boundaries',
});

export const showPopulationState = selector({
  key: 'showPopulationState',
  get: ({ get }) => get(regionDataState) === 'population',
});
