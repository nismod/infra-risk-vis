import React from 'react';
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormControl from '@material-ui/core/FormControl';
import FormLabel from '@material-ui/core/FormLabel';

const BackgroundControl = (props) => (
  <FormControl component="fieldset">
    <FormLabel component="legend">Map Background</FormLabel>
    <RadioGroup
      aria-label="Map background"
      name="map_background"
      value={props.background}
      onChange={props.onBackgroundChange}
    >
      <FormControlLabel value="light" control={<Radio />} label="Map" />
      <FormControlLabel value="satellite" control={<Radio />} label="Satellite" />
    </RadioGroup>
  </FormControl>
);

export default BackgroundControl;
