import { FC, useEffect, useState } from 'react';
import { Box, Checkbox, FormControl, FormControlLabel, FormLabel, makeStyles, Typography } from '@material-ui/core';

import { CustomNumberSlider } from './CustomSlider';
import { getHazardId, LAYERS } from '../config/layers';

interface HazardsControlProps {
  layerVisibility: Record<string, boolean>;
  onLayerVisibilityUpdate: (visibilityUpdate: Record<string, boolean>) => void;
}

const fluvialLayers = Object.keys(LAYERS).filter((l) => l.startsWith('fluvial'));

const fluvialMarks = [
  { value: 20, label: '20' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
  { value: 200, label: '200' },
  { value: 500, label: '500' },
  { value: 1500, label: '1500' },
];

const coastalLayers = Object.keys(LAYERS).filter((l) => l.startsWith('coastal'));
const coastalMarks = [
  { value: 1, label: '1' },
  { value: 2, label: '2' },
  { value: 5, label: '5' },
  { value: 10, label: '10' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
];

const surfaceLayers = Object.keys(LAYERS).filter((l) => l.startsWith('surface'));
const surfaceMarks = [
  { value: 20, label: '20' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
  { value: 200, label: '200' },
  { value: 500, label: '500' },
  { value: 1500, label: '1500' },
];

const useStyles = makeStyles({
  sliderForm: {
    width: '100%',
  },
});

export const HazardsControl: FC<HazardsControlProps> = ({ layerVisibility, onLayerVisibilityUpdate }) => {
  const classes = useStyles();

  const [showFluvial, setShowFluvial] = useState(false);
  const [fluvialReturnPeriod, setFluvialReturnPeriod] = useState(20);

  // TODO modify to avoid using useEffect with incomplete dependencies
  useEffect(() => {
    const layerToShow = getHazardId({
      hazardType: 'fluvial',
      returnPeriod: fluvialReturnPeriod,
      rcp: 'baseline',
      epoch: 2010,
      confidence: 'None',
    });
    onLayerVisibilityUpdate(Object.fromEntries(fluvialLayers.map((ln) => [ln, showFluvial && ln === layerToShow])));
  }, [fluvialReturnPeriod, showFluvial]);

  const [showCoastal, setShowCoastal] = useState(false);
  const [coastalReturnPeriod, setCoastalReturnPeriod] = useState(1);

  // TODO modify to avoid using useEffect with incomplete dependencies
  useEffect(() => {
    const layerToShow = getHazardId({
      hazardType: 'coastal',
      returnPeriod: coastalReturnPeriod,
      rcp: '4x5',
      epoch: 2050,
      confidence: 'None',
    });
    // const layerToShow = `hazard_coastal_${coastalReturnPeriod}`;
    onLayerVisibilityUpdate(Object.fromEntries(coastalLayers.map((ln) => [ln, showCoastal && ln === layerToShow])));
  }, [showCoastal, coastalReturnPeriod]);

  const [showSurface, setShowSurface] = useState(false);
  const [surfaceReturnPeriod, setSurfaceReturnPeriod] = useState(20);

  // TODO modify to avoid using useEffect with incomplete dependencies
  useEffect(() => {
    const layerToShow = getHazardId({
      hazardType: 'surface',
      returnPeriod: surfaceReturnPeriod,
      rcp: 'baseline',
      epoch: 2010,
      confidence: 'None',
    });
    onLayerVisibilityUpdate(Object.fromEntries(surfaceLayers.map((ln) => [ln, showSurface && ln === layerToShow])));
  }, [showSurface, surfaceReturnPeriod]);

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
            label="River Flooding"
          />
        </FormControl>
        <FormControl disabled={!showFluvial} component="fieldset" className={classes.sliderForm}>
          <FormLabel component="legend">Return Period</FormLabel>
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
            label="Coastal Flooding"
          />
        </FormControl>
        <FormControl disabled={!showCoastal} component="fieldset" className={classes.sliderForm}>
          <FormLabel component="legend">Return Period</FormLabel>
          <CustomNumberSlider
            marks={coastalMarks}
            value={coastalReturnPeriod}
            onChange={setCoastalReturnPeriod}
            disabled={!showCoastal}
          />
        </FormControl>
      </Box>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={showSurface}
                onChange={(e) => setShowSurface(e.currentTarget.checked)}
              />
            }
            label="Surface Flooding"
          />
        </FormControl>
        <FormControl disabled={!showSurface} component="fieldset" className={classes.sliderForm}>
          <FormLabel component="legend">Return Period</FormLabel>
          <CustomNumberSlider
            marks={surfaceMarks}
            value={surfaceReturnPeriod}
            onChange={setSurfaceReturnPeriod}
            disabled={!showSurface}
          />
        </FormControl>
      </Box>
    </Box>
  );
};
