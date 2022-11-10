import { FormControl, InputLabel, MenuItem, Select } from '@mui/material';
import { selector, useRecoilState, useRecoilValue } from 'recoil';

import { titleCase, unique } from '@/lib/helpers';
import { makeSelectState } from '@/lib/recoil/make-state/make-select-state';

import { HAZARDS_METADATA } from '@/config/hazards/metadata';

import { damagesDataState } from './ExpectedDamagesSection';

export const hazardsState = selector({
  key: 'DamagesSection/hazardsState',
  get: ({ get }) => unique(get(damagesDataState).map((d) => d.hazard)),
});
export const selectedHazardState = makeSelectState('DamagesSection/selectedHazard', hazardsState);

export const HazardSelect = () => {
  const hazards = useRecoilValue(hazardsState);
  const [selectedHazard, setSelectedHazard] = useRecoilState(selectedHazardState);

  return hazards.length ? (
    <FormControl fullWidth sx={{ my: 2 }} disabled={hazards.length === 1}>
      <InputLabel>Hazard</InputLabel>
      <Select label="Hazard" value={selectedHazard ?? ''} onChange={(e) => setSelectedHazard(e.target.value as string)}>
        {hazards.map((h) => (
          <MenuItem key={h} value={h}>
            {HAZARDS_METADATA[h]?.label ?? titleCase(h)}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  ) : null;
};

export const epochsState = selector({
  key: 'DamagesSection/epochsState',
  get: ({ get }) => unique(get(damagesDataState).map((d) => d.epoch)).sort(),
});
export const selectedEpochState = makeSelectState('DamagesSection/selectedEpoch', epochsState);

export const EpochSelect = () => {
  const epochs = useRecoilValue(epochsState);
  const [selectedEpoch, setSelectedEpoch] = useRecoilState(selectedEpochState);

  return epochs.length ? (
    <FormControl fullWidth disabled={epochs.length === 1}>
      <InputLabel>Epoch</InputLabel>
      <Select label="Epoch" value={selectedEpoch ?? ''} onChange={(e) => setSelectedEpoch(e.target.value as string)}>
        {epochs.map((h) => (
          <MenuItem key={h} value={h}>
            {titleCase(h)}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  ) : null;
};
