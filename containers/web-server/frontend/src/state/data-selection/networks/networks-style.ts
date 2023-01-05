import { selector } from 'recoil';

import { sidebarPathVisibilityState } from '@/sidebar/SidebarContent';
import { viewState } from '@/state/view';

export const networksStyleState = selector<string>({
  key: 'networksStyleState',
  get: ({ get }) =>
    get(viewState) === 'risk' && get(sidebarPathVisibilityState('risk/infrastructure'))
      ? 'damages'
      : 'type',
});
