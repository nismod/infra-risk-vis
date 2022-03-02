import { selector } from 'recoil';
import { showDirectDamagesState } from './damage-mapping/damage-map';
import { damageMapStyleParamsState } from './damage-mapping/damage-style-params';

export const styleParamsState = selector({
  key: 'styleParamsState',
  get: ({ get }) => {
    if (!get(showDirectDamagesState)) return {};

    return get(damageMapStyleParamsState);
  },
});
