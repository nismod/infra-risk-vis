import { NETWORK_LAYER_NAMES } from 'config/networks/metadata';
import { atom } from 'recoil';

export const networkSelectionState = atom({
  key: 'networkSelectionState',
  default: Object.fromEntries(NETWORK_LAYER_NAMES.map((nl) => [nl, false])),
});
