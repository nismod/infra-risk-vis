import { atom, selector, selectorFamily } from 'recoil';
import { networksStyleState } from 'state/networks/networks-style';

export const showDirectDamagesState = selector({
  key: 'showDamagesState',
  get: ({ get }) => get(networksStyleState) === 'direct-damages',
});

export const selectedDamageSourceState = atom({
  key: 'damageSourceSelectionImpl',
  default: 'total-damages',
});

export const damageSourceSelectionState = selectorFamily<boolean, string>({
  key: 'damageSourceSelectionState',
  get:
    (damageSourceId) =>
    ({ get }) =>
      get(selectedDamageSourceState) === damageSourceId,
  set:
    (damageSourceId) =>
    ({ set }, newValue) => {
      if (newValue !== false) {
        set(selectedDamageSourceState, damageSourceId);
      }
    },
});
