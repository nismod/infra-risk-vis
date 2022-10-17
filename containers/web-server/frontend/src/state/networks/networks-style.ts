import { selector } from 'recoil';

import { sectionStyleValueState } from '@/state/sections';

export const networksStyleState = selector({
  key: 'networksStyleState',
  get: ({ get }) => get(sectionStyleValueState('assets')),
});
