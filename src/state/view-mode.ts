import { atom } from 'recoil';

export const viewModeState = atom<'input' | 'direct-damages'>({
  key: 'viewModeState',
  default: 'input',
});
