import { NETWORKS_DEFAULT_STYLE } from 'config/networks/styles';
import { atom } from 'recoil';

export const networksStyleState = atom({
  key: 'networksStyleState',
  default: NETWORKS_DEFAULT_STYLE,
});
