import React, { FC } from 'react';
import { FormControl, FormControlLabel, FormLabel, Radio, RadioGroup } from '@mui/material';
import { BackgroundName } from '../config/backgrounds';

interface BackgroundControlProps {
  background: BackgroundName;
  onBackgroundChange: (background: BackgroundName) => void;
}

export const BackgroundControl: FC<BackgroundControlProps> = ({ background, onBackgroundChange }) => (
  <FormControl component="fieldset">
    <FormLabel component="legend">Map Background</FormLabel>
    <RadioGroup
      aria-label="Map background"
      name="map_background"
      value={background}
      onChange={(e) => onBackgroundChange(e.target.value as BackgroundName)}
    >
      <FormControlLabel value="light" control={<Radio />} label="Map" />
      <FormControlLabel value="satellite" control={<Radio />} label="Satellite" />
    </RadioGroup>
  </FormControl>
);
