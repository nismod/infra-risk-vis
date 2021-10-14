import { FC, useEffect, useState } from 'react';
import {
  Box,
  Checkbox,
  FormControl,
  FormControlLabel,
  FormLabel,
  MenuItem,
  Select,
  Typography,
} from '@material-ui/core';

import { CustomNumberSlider } from './CustomSlider';
import { getHazardId, LAYERS } from '../config/layers';

interface HazardsControlProps {
  layerVisibility: Record<string, boolean>;
  onLayerVisibilityUpdate: (visibilityUpdate: Record<string, boolean>) => void;
}

const rcpLookup = {
  baseline: 'Baseline',
  '2x6': '2.6',
  '4x5': '4.5',
  '8x5': '8.5',
};

function rcpLabel(value) {
  return rcpLookup[value];
}

function epochLabel(value) {
  if (value === 2010) return 'Present';
  return value;
}

// const fluvialLayers = Object.keys(LAYERS).filter((l) => l.startsWith('fluvial'));

// const fluvialMarks = [
//   { value: 20, label: '20' },
//   { value: 50, label: '50' },
//   { value: 100, label: '100' },
//   { value: 200, label: '200' },
//   { value: 500, label: '500' },
//   { value: 1500, label: '1500' },
// ];

// const coastalLayers = Object.keys(LAYERS).filter((l) => l.startsWith('coastal'));
// const coastalMarks = [
//   { value: 1, label: '1' },
//   { value: 2, label: '2' },
//   { value: 5, label: '5' },
//   { value: 10, label: '10' },
//   { value: 50, label: '50' },
//   { value: 100, label: '100' },
// ];

// const surfaceLayers = Object.keys(LAYERS).filter((l) => l.startsWith('surface'));
// const surfaceMarks = [
//   { value: 20, label: '20' },
//   { value: 50, label: '50' },
//   { value: 100, label: '100' },
//   { value: 200, label: '200' },
//   { value: 500, label: '500' },
//   { value: 1500, label: '1500' },
// ];

export const HazardsControl = ({
  hazardShow,
  hazardOptions,
  hazardParams,
  onSingleHazardShow,
  onSingleHazardParam,
}) => {
  // const [showFluvial, setShowFluvial] = useState(false);
  // const [fluvialReturnPeriod, setFluvialReturnPeriod] = useState(20);

  // // TODO modify to avoid using useEffect with incomplete dependencies
  // useEffect(() => {
  //   const layerToShow = getHazardId({
  //     hazardType: 'fluvial',
  //     returnPeriod: fluvialReturnPeriod,
  //     rcp: 'baseline',
  //     epoch: 2010,
  //     confidence: 'None',
  //   });
  //   onLayerVisibilityUpdate(Object.fromEntries(fluvialLayers.map((ln) => [ln, showFluvial && ln === layerToShow])));
  // }, [fluvialReturnPeriod, showFluvial]);

  // const [showCoastal, setShowCoastal] = useState(false);
  // const [coastalReturnPeriod, setCoastalReturnPeriod] = useState(1);

  // // TODO modify to avoid using useEffect with incomplete dependencies
  // useEffect(() => {
  //   const layerToShow = getHazardId({
  //     hazardType: 'coastal',
  //     returnPeriod: coastalReturnPeriod,
  //     rcp: '4x5',
  //     epoch: 2050,
  //     confidence: 'None',
  //   });
  //   // const layerToShow = `hazard_coastal_${coastalReturnPeriod}`;
  //   onLayerVisibilityUpdate(Object.fromEntries(coastalLayers.map((ln) => [ln, showCoastal && ln === layerToShow])));
  // }, [showCoastal, coastalReturnPeriod]);

  // const [showSurface, setShowSurface] = useState(false);
  // const [surfaceReturnPeriod, setSurfaceReturnPeriod] = useState(20);

  // // TODO modify to avoid using useEffect with incomplete dependencies
  // useEffect(() => {
  //   const layerToShow = getHazardId({
  //     hazardType: 'surface',
  //     returnPeriod: surfaceReturnPeriod,
  //     rcp: 'baseline',
  //     epoch: 2010,
  //     confidence: 'None',
  //   });
  //   onLayerVisibilityUpdate(Object.fromEntries(surfaceLayers.map((ln) => [ln, showSurface && ln === layerToShow])));
  // }, [showSurface, surfaceReturnPeriod]);

  return (
    <Box mb={2}>
      <Typography variant="h6">Hazard Layers</Typography>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={hazardShow.fluvial}
                onChange={(e) => onSingleHazardShow('fluvial', e.currentTarget.checked)}
              />
            }
            label="River Flooding"
          />
        </FormControl>
        <FormControl disabled={!hazardShow.fluvial} component="fieldset" fullWidth>
          <FormLabel component="legend">Return Period</FormLabel>
          <CustomNumberSlider
            marks={hazardOptions.fluvial.returnPeriod}
            value={hazardParams.fluvial.returnPeriod}
            onChange={(v) => onSingleHazardParam('fluvial', 'returnPeriod', v)}
            disabled={!hazardShow.fluvial}
          />
        </FormControl>
      </Box>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={hazardShow.coastal}
                onChange={(e) => onSingleHazardShow('coastal', e.currentTarget.checked)}
              />
            }
            label="Coastal Flooding"
          />
        </FormControl>
        <Box mb={2}>
          <FormControl disabled={!hazardShow.coastal} component="fieldset" fullWidth>
            <FormLabel component="legend">Return Period</FormLabel>
            <CustomNumberSlider
              marks={hazardOptions.coastal.returnPeriod}
              value={hazardParams.coastal.returnPeriod}
              onChange={(v) => onSingleHazardParam('coastal', 'returnPeriod', v)}
              disabled={!hazardShow.coastal}
            />
          </FormControl>
        </Box>
        <Box mb={2}>
          <FormControl disabled={!hazardShow.coastal} component="fieldset" fullWidth>
            <FormLabel component="legend">RCP</FormLabel>
            <Select
              value={hazardParams.coastal.rcp}
              onChange={(e) => onSingleHazardParam('coastal', 'rcp', e.target.value)}
            >
              {hazardOptions.coastal.rcp.map((rcp) => (
                <MenuItem value={rcp}>{rcpLookup[rcp]}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Box>
        <Box mb={2}>
          <FormControl disabled={!hazardShow.coastal} component="fieldset" fullWidth>
            <FormLabel component="legend">Epoch</FormLabel>
            <Select
              value={hazardParams.coastal.epoch}
              onChange={(e) => onSingleHazardParam('coastal', 'epoch', e.target.value)}
            >
              {hazardOptions.coastal.epoch.map((epoch) => (
                <MenuItem value={epoch}>{epochLabel(epoch)}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Box>
      </Box>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={hazardShow.surface}
                onChange={(e) => onSingleHazardShow('surface', e.currentTarget.checked)}
              />
            }
            label="Surface Flooding"
          />
        </FormControl>
        <FormControl disabled={!hazardShow.surface} component="fieldset" fullWidth>
          <FormLabel component="legend">Return Period</FormLabel>
          <CustomNumberSlider
            marks={hazardOptions.surface.returnPeriod}
            value={hazardParams.surface.returnPeriod}
            onChange={(v) => onSingleHazardParam('surface', 'returnPeriod', v)}
            disabled={!hazardShow.surface}
          />
        </FormControl>
      </Box>
      <Box mb={1}>
        <FormControl>
          <FormControlLabel
            control={
              <Checkbox
                color="primary"
                checked={hazardShow.cyclone}
                onChange={(e) => onSingleHazardShow('cyclone', e.currentTarget.checked)}
              />
            }
            label="Cyclones"
          />
        </FormControl>
        <FormControl disabled={!hazardShow.surface} component="fieldset" fullWidth>
          <FormLabel component="legend">Return Period</FormLabel>
          <Select
            value={hazardParams.cyclone.returnPeriod}
            onChange={(e) => onSingleHazardParam('cyclone', 'returnPeriod', e.target.value)}
            disabled={!hazardShow.cyclone}
          >
            {hazardOptions.cyclone.returnPeriod.map((rp) => (
              <MenuItem value={rp}>{rp}</MenuItem>
            ))}
          </Select>
        </FormControl>
      </Box>
    </Box>
  );
};
