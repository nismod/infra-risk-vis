import {
  Box,
  Checkbox,
  FormControl,
  FormControlLabel,
  FormLabel,
  InputLabel,
  MenuItem,
  Select,
  Typography,
} from '@material-ui/core';

import { CustomNumberSlider } from './CustomSlider';

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

const HazardSection = ({ children }) => <Box mb={4}>{children}</Box>;
const InputSection = ({ children }) => <Box mb={2}>{children}</Box>;

export const HazardsControl = ({
  hazardShow,
  hazardOptions,
  hazardParams,
  onSingleHazardShow,
  onSingleHazardParam,
}) => {
  return (
    <Box mb={1}>
      <Typography variant="h6">Hazard Layers</Typography>
      <HazardSection>
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
        <FormControl disabled={!hazardShow.fluvial} fullWidth>
          <FormLabel>Return Period</FormLabel>
          <CustomNumberSlider
            marks={hazardOptions.fluvial.returnPeriod}
            value={hazardParams.fluvial.returnPeriod}
            onChange={(v) => onSingleHazardParam('fluvial', 'returnPeriod', v)}
            disabled={!hazardShow.fluvial}
          />
        </FormControl>
      </HazardSection>
      <HazardSection>
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
        <FormControl disabled={!hazardShow.surface} fullWidth>
          <FormLabel>Return Period</FormLabel>
          <CustomNumberSlider
            marks={hazardOptions.surface.returnPeriod}
            value={hazardParams.surface.returnPeriod}
            onChange={(v) => onSingleHazardParam('surface', 'returnPeriod', v)}
            disabled={!hazardShow.surface}
          />
        </FormControl>
      </HazardSection>
      <HazardSection>
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
        </Box>
        <InputSection>
          <FormControl disabled={!hazardShow.coastal} component="fieldset" fullWidth>
            <FormLabel>Return Period</FormLabel>
            <CustomNumberSlider
              marks={hazardOptions.coastal.returnPeriod}
              value={hazardParams.coastal.returnPeriod}
              onChange={(v) => onSingleHazardParam('coastal', 'returnPeriod', v)}
              disabled={!hazardShow.coastal}
            />
          </FormControl>
        </InputSection>
        <InputSection>
          <FormControl disabled={!hazardShow.coastal} variant="standard">
            <InputLabel>Epoch</InputLabel>
            <Select
              value={hazardParams.coastal.epoch}
              onChange={(e) => onSingleHazardParam('coastal', 'epoch', e.target.value)}
            >
              {hazardOptions.coastal.epoch.map((epoch) => (
                <MenuItem value={epoch}>{epochLabel(epoch)}</MenuItem>
              ))}
            </Select>
          </FormControl>
          <FormControl disabled={!hazardShow.coastal}>
            <InputLabel>RCP</InputLabel>
            <Select
              value={hazardParams.coastal.rcp}
              onChange={(e) => onSingleHazardParam('coastal', 'rcp', e.target.value)}
            >
              {hazardOptions.coastal.rcp.map((rcp) => (
                <MenuItem value={rcp}>{rcpLabel(rcp)}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </InputSection>
      </HazardSection>
      <HazardSection>
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
        <InputSection>
          <FormControl disabled={!hazardShow.cyclone} component="fieldset" fullWidth>
            <FormLabel>Return Period</FormLabel>
            <CustomNumberSlider
              marks={hazardOptions.cyclone.returnPeriod}
              value={hazardParams.cyclone.returnPeriod}
              onChange={(v) => onSingleHazardParam('cyclone', 'returnPeriod', v)}
              disabled={!hazardShow.cyclone}
              valueLabelDisplay="auto"
              showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
            />
          </FormControl>
        </InputSection>
        <InputSection>
          <FormControl disabled={!hazardShow.cyclone} variant="standard">
            <InputLabel>Epoch</InputLabel>
            <Select
              value={hazardParams.cyclone.epoch}
              onChange={(e) => onSingleHazardParam('cyclone', 'epoch', e.target.value)}
            >
              {hazardOptions.cyclone.epoch.map((epoch) => (
                <MenuItem value={epoch}>{epochLabel(epoch)}</MenuItem>
              ))}
            </Select>
          </FormControl>
          <FormControl disabled={!hazardShow.cyclone}>
            <InputLabel>RCP</InputLabel>
            <Select
              value={hazardParams.cyclone.rcp}
              onChange={(e) => onSingleHazardParam('cyclone', 'rcp', e.target.value)}
            >
              {hazardOptions.cyclone.rcp.map((rcp) => (
                <MenuItem value={rcp}>{rcpLabel(rcp)}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </InputSection>
      </HazardSection>
    </Box>
  );
};
