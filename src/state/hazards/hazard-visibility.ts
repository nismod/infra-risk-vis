import { HAZARD_LAYER_NAMES } from 'config/hazards/metadata';
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
      return _.fromPairs(HAZARD_LAYER_NAMES.map((group) => [group, get(hazardSelectionState(group))]));
    }
  },
});
