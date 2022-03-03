import { atomFamily } from 'recoil';

export const showLayerState = atomFamily<boolean, string>({
  key: 'showLayerState',
  default: true,
});
