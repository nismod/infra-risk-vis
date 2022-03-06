import { RegionLevel } from 'config/regions/metadata';
import { REGION_DEFAULT_STYLE } from 'config/regions/styles';
import { atom, selector } from 'recoil';

export const regionLevelState = atom<RegionLevel>({
  key: 'regionLevelState',
  default: 'parish',
});

export const regionsStyleState = atom({
  key: 'regionDataState',
  default: REGION_DEFAULT_STYLE,
});

export const showPopulationState = selector({
  key: 'showPopulationState',
  get: ({ get }) => get(regionsStyleState) === 'population',
});
