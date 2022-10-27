import { atom, selector } from 'recoil';

import { RegionLevel } from '@/config/_old/regions/metadata';

import { sectionStyleValueState } from '../../sections';

export const regionLevelState = atom<RegionLevel>({
  key: 'regionLevelState',
  default: 'parish',
});

export const regionsStyleState = selector({
  key: 'regionsStyleState',
  get: ({ get }) => get(sectionStyleValueState('regions')),
});

export const showPopulationState = selector({
  key: 'showPopulationState',
  get: ({ get }) => get(regionsStyleState) === 'population',
});
