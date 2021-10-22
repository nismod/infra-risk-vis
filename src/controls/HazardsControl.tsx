import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Box,
  Checkbox,
  FormControl,
  FormControlLabel,
  FormLabel,
  Grid,
  InputLabel,
  MenuItem,
  Select,
  Switch,
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

const HazardSection = ({ show, onShow, label, children }) => (
  <Accordion expanded={show} onChange={onShow}>
    <AccordionSummary>
      <Typography>{label}</Typography>
    </AccordionSummary>
    <AccordionDetails style={{ display: 'block' }}>{children}</AccordionDetails>
  </Accordion>
);
const InputSection = ({ children }) => (
  <Box mb={2} flexGrow={1} width="100%">
    {children}
  </Box>
);

export const HazardsControl = ({
  hazardShow,
  hazardOptions,
  hazardParams,
  onSingleHazardShow,
  onSingleHazardParam,
  showDamages,
  showDamageRaster,
  onShowDamageRaster,
}) => {
  const handleChange = (hazardType) => (e, isExpanded) => {
    onSingleHazardShow(hazardType, isExpanded);
  };

  return (
    <Box mb={1}>
      <Box mt={2} mb={1}>
        <Grid container justifyContent="space-between">
          <Grid item>
            <Typography variant="h6">Hazards</Typography>
          </Grid>
          {showDamages && (
            <Grid item>
              <FormControlLabel
                control={
                  <Switch checked={showDamageRaster} onChange={(e, checked) => onShowDamageRaster(checked)}></Switch>
                }
                label="Show hazard raster"
                labelPlacement="start"
              />
            </Grid>
          )}
        </Grid>
      </Box>
      <HazardSection show={hazardShow.fluvial} onShow={handleChange('fluvial')} label="River Flooding">
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
      <HazardSection show={hazardShow.surface} onShow={handleChange('surface')} label="Surface Flooding">
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
      <HazardSection show={hazardShow.coastal} onShow={handleChange('coastal')} label="Coastal Flooding">
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
          <FormControl disabled={!hazardShow.coastal} variant="standard" style={{ width: '50%' }}>
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
          <FormControl disabled={!hazardShow.coastal} style={{ width: '50%' }}>
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
      <HazardSection show={hazardShow.cyclone} onShow={handleChange('cyclone')} label="Cyclones">
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
          <FormControl disabled={!hazardShow.cyclone} variant="standard" style={{ width: '50%' }}>
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
          <FormControl disabled={!hazardShow.cyclone} style={{ width: '50%' }}>
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
