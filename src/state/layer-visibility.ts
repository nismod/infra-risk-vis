import { selector } from 'recoil';
import { hazardVisibilityState } from './hazards/hazard-visibility';
import { networkSelectionState } from './network-selection';

export const layerVisibilityState = selector({
  key: 'layerVisibilityState',
  get: ({ get }) => Object.assign({}, get(networkSelectionState), get(hazardVisibilityState)),
});
