import { atomFamily } from 'recoil';

export const hazardSelectionState = atomFamily({
  key: 'hazardSelectionState',
  default: false,
});
