import { hazardTypes } from 'config/data/hazards';
import _ from 'lodash';
import { selector } from 'recoil';
import {
  selectedDamageSourceState,
  showDamageRasterState,
  showDirectDamagesState,
} from 'state/damage-mapping/damage-map';
import { hazardSelectionState } from './hazard-selection';

export const hazardVisibilityState = selector({
  key: 'hazardVisibilityState',
  get: ({ get }) => {
    if (get(showDirectDamagesState)) {
      const selectedDamageSource = get(selectedDamageSourceState);
      if (selectedDamageSource === 'total-damages' || !get(showDamageRasterState)) {
        return {};
      } else {
        return {
          [selectedDamageSource]: true,
        };
      }
    } else {
      return _.fromPairs(hazardTypes.map((group) => [group, get(hazardSelectionState(group))]));
    }
  },
});
