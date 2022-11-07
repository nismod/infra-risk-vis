import { selector } from 'recoil';

import { sectionStyleValueState } from '@/state/sections';
import { viewState } from '@/state/view';

export const networksStyleState = selector<string>({
  key: 'networksStyleState',
  get: ({ get }) => (get(viewState) === 'risk' ? 'damages' : 'type'), //get(sectionStyleValueState('assets')),
});
