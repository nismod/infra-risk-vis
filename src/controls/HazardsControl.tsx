import { FC, useEffect, useState } from 'react';
import { Box, Checkbox, FormControl, FormControlLabel, FormLabel, makeStyles, Typography } from '@material-ui/core';

import { CustomNumberSlider } from './CustomSlider';
import { LAYERS } from '../config/layers';

interface HazardsControlProps {
  layerVisibility: Record<string, boolean>;
  onLayerVisibilityUpdate: (visibilityUpdate: Record<string, boolean>) => void;
}

const fluvialLayers = Object.keys(LAYERS).filter((l) => l.startsWith('hazard_fluvial'));

const fluvialMarks = [
  { value: 20, label: '20' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
  { value: 200, label: '200' },
  { value: 500, label: '500' },
  { value: 1500, label: '1500' },
];

const coastalLayers = Object.keys(LAYERS).filter((l) => l.startsWith('hazard_coastal'));
const coastalMarks = [
  { value: 1, label: '1' },
  { value: 2, label: '2' },
  { value: 5, label: '5' },
  { value: 10, label: '10' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
];

const useStyles = makeStyles({
  sliderForm: {
    width: '100%',
  },
});

export const HazardsControl: FC<HazardsControlProps> = ({ layerVisibility, onLayerVisibilityUpdate }) => {
  const classes = useStyles();

  const [showFluvial, setShowFluvial] = useState(false);
  const [fluvialReturnPeriod, setFluvialReturnPeriod] = useState(10);

  const [showCoastal, setShowCoastal] = useState(false);
  const [coastalReturnPeriod, setCoastalReturnPeriod] = useState(1);

  // TODO modify to avoid using useEffect with incomplete dependencies
  useEffect(() => {
    const layerToShow = `hazard_fluvial_${fluvialReturnPeriod}`;
    onLayerVisibilityUpdate(Object.fromEntries(fluvialLayers.map((ln) => [ln, showFluvial && ln === layerToShow])));
  }, [fluvialReturnPeriod, showFluvial]);

  // TODO modify to avoid using useEffect with incomplete dependencies
  useEffect(() => {
    const layerToShow = `hazard_coastal_${coastalReturnPeriod}`;
    onLayerVisibilityUpdate(Object.fromEntries(coastalLayers.map((ln) => [ln, showCoastal && ln === layerToShow])));
  }, [showCoastal, coastalReturnPeriod]);

  return (
    <Box mb={2}>
      <Typography variant="h6">Hazard Layers</Typography>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={showFluvial}
                onChange={(e) => setShowFluvial(e.currentTarget.checked)}
              />
            }
            label="Show River Flooding"
          />
        </FormControl>
        <FormControl disabled={!showFluvial} component="fieldset" className={classes.sliderForm}>
          <FormLabel component="legend">River Flooding Return Period</FormLabel>
          <CustomNumberSlider
            marks={fluvialMarks}
            value={fluvialReturnPeriod}
            onChange={setFluvialReturnPeriod}
            disabled={!showFluvial}
          />
        </FormControl>
      </Box>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={showCoastal}
                onChange={(e) => setShowCoastal(e.currentTarget.checked)}
              />
            }
            label="Show Coastal Flooding"
          />
        </FormControl>
        <FormControl disabled={!showCoastal} component="fieldset" className={classes.sliderForm}>
          <FormLabel component="legend">Coastal Flooding Return Period</FormLabel>
          <CustomNumberSlider
            marks={coastalMarks}
            value={coastalReturnPeriod}
            onChange={setCoastalReturnPeriod}
            disabled={!showCoastal}
          />
        </FormControl>
      </Box>
    </Box>
  );
};
