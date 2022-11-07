import { selector } from 'recoil';

import { sidebarVisibilityToggleState } from '@/sidebar/SidebarContent';
import { viewState } from '@/state/view';

export const networksStyleState = selector<string>({
  key: 'networksStyleState',
  get: ({ get }) =>
    get(viewState) === 'risk' && get(sidebarVisibilityToggleState('risk/infrastructure')) ? 'damages' : 'type',
});
