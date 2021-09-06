import React, { FC } from 'react';
import { FormControl, FormControlLabel, FormLabel, Radio, RadioGroup } from '@material-ui/core';

interface BackgroundControlProps {
  background: string;
  onBackgroundChange: (background: string) => void;
}

export const BackgroundControl: FC<BackgroundControlProps> = ({ background, onBackgroundChange }) => (
  <FormControl component="fieldset">
    <FormLabel component="legend">Map Background</FormLabel>
    <RadioGroup
      aria-label="Map background"
      name="map_background"
      value={background}
      onChange={(e) => onBackgroundChange(e.target.value)}
    >
      <FormControlLabel value="light" control={<Radio />} label="Map" />
      <FormControlLabel value="satellite" control={<Radio />} label="Satellite" />
    </RadioGroup>
  </FormControl>
);
