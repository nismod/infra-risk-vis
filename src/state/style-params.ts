import { selector } from 'recoil';
import { showDamagesState } from './damage-mapping/damage-map';
import { damageMapStyleParamsState } from './damage-mapping/damage-style-params';

export const styleParamsState = selector({
  key: 'styleParamsState',
  get: ({ get }) => {
    if (!get(showDamagesState)) return {};

    return get(damageMapStyleParamsState);
  },
});
