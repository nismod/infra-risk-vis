import { selector } from 'recoil';
import { sectionStyleValueState, sectionVisibilityState } from './sections';

export const showAdaptationsTableState = selector<boolean>({
  key: 'showAdaptationsTable',
  get: ({ get }) => get(sectionVisibilityState('assets')) && get(sectionStyleValueState('assets')) === 'adaptation',
});
