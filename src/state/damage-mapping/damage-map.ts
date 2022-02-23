import { atom, selector, selectorFamily } from 'recoil';
import { viewModeState } from '../view-mode';

export const showDirectDamagesState = selector({
  key: 'showDamagesState',
  get: ({ get }) => get(viewModeState) === 'direct-damages',
});

export const showDamageRasterState = atom({
  key: 'showDamageRasterState',
  default: true,
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
    ({ get, set }, newValue) => {
      const currentDamageSourceId = get(selectedDamageSourceState);
      if (newValue === false) {
        if (damageSourceId === currentDamageSourceId) {
          set(selectedDamageSourceState, null);
        }
      } else {
        set(selectedDamageSourceState, damageSourceId);
      }
    },
});
