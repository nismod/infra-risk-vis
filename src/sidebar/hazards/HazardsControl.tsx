import { Box, FormControl, FormLabel, Grid, MenuItem, Select } from '@mui/material';
import { Children } from 'react';
import { useRecoilValue } from 'recoil';

import { ToggleSection, ToggleSectionGroup } from 'lib/controls/accordion-toggle/ToggleSection';
import { CustomNumberSlider } from 'lib/controls/CustomSlider';

import { damageSourceSelectionState, showDirectDamagesState } from 'state/damage-mapping/damage-map';
import { dataParamOptionsState, dataParamState, useUpdateDataParam } from '../../state/data-params';
import { hazardSelectionState } from '../../state/hazards/hazard-selection';

const DataParam = ({ group, id, children }) => {
  const value = useRecoilValue(dataParamState({ group, param: id }));
  const updateValue = useUpdateDataParam(group, id);
  const options = useRecoilValue(dataParamOptionsState({ group, param: id }));

  return typeof children === 'function' ? children({ value: value, onChange: updateValue, options }) : children;
};

const InputSection = ({ children }) => (
  <Box mb={2} flexGrow={1} width="100%">
    {children}
  </Box>
);

const ReturnPeriodControl = ({ group, ...otherProps }) => {
  return (
    <FormControl fullWidth>
      <FormLabel>Return Period</FormLabel>
      <DataParam group={group} id="returnPeriod">
        {({ value, onChange, options }) => (
          <CustomNumberSlider marks={options} value={value} onChange={onChange} {...otherProps} />
        )}
      </DataParam>
    </FormControl>
  );
};

function epochLabel(value) {
  if (value === 2010) return 'Present';
  return value;
}

const EpochControl = ({ group }) => {
  return (
    <FormControl fullWidth>
      <FormLabel>Epoch</FormLabel>
      <DataParam group={group} id="epoch">
        {({ value, onChange, options }) => (
          <Select variant="standard" value={value} onChange={(e) => onChange(e.target.value)} fullWidth>
            {options.map((epoch) => (
              <MenuItem key={epoch} value={epoch}>
                {epochLabel(epoch)}
              </MenuItem>
            ))}
          </Select>
        )}
      </DataParam>
    </FormControl>
  );
};

function rcpLabel(value) {
  return value === 'baseline' ? 'Baseline' : value;
}

const RCPControl = ({ group }) => {
  return (
    <FormControl fullWidth>
      <FormLabel>RCP</FormLabel>
      <DataParam group={group} id="rcp">
        {({ value, onChange, options }) => (
          <Select variant="standard" value={value} onChange={(e) => onChange(e.target.value)} fullWidth>
            {options.map((rcp) => (
              <MenuItem key={rcp} value={rcp}>
                {rcpLabel(rcp)}
              </MenuItem>
            ))}
          </Select>
        )}
      </DataParam>
    </FormControl>
  );
};

const ParamRow = ({ children }) => (
  <Grid container spacing={1}>
    {Children.map(children, (child) => (
      <Grid item xs>
        {child}
      </Grid>
    ))}
  </Grid>
);

export const HazardsControl = () => {
  const showDamages = useRecoilValue(showDirectDamagesState);
  const forceSingle = showDamages;

  const selectionState = showDamages ? damageSourceSelectionState : hazardSelectionState;

  return (
    <Box mb={1}>
      <ToggleSectionGroup toggleState={selectionState}>
        {showDamages && (
          <ToggleSection id="total-damages" label="Total Damages" forceSingle={forceSingle}>
            <InputSection>
              <ParamRow>
                <EpochControl group="total-damages" />
                <RCPControl group="total-damages" />
              </ParamRow>
            </InputSection>
          </ToggleSection>
        )}

        <ToggleSection id="fluvial" label="River Flooding" forceSingle={forceSingle}>
          <ReturnPeriodControl group="fluvial" param="returnPeriod" />
        </ToggleSection>

        <ToggleSection id="surface" label="Surface Flooding" forceSingle={forceSingle}>
          <ReturnPeriodControl group="surface" param="returnPeriod" />
        </ToggleSection>

        <ToggleSection id="coastal" label="Coastal Flooding" forceSingle={forceSingle}>
          <InputSection>
            <ReturnPeriodControl group="coastal" param="returnPeriod" />
          </InputSection>
          <InputSection>
            <ParamRow>
              <EpochControl group="coastal" />
              <RCPControl group="coastal" />
            </ParamRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="cyclone" label="Cyclones" forceSingle={forceSingle}>
          <InputSection>
            <ReturnPeriodControl
              group="cyclone"
              valueLabelDisplay="auto"
              showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
            />
          </InputSection>
          <InputSection>
            <ParamRow>
              <EpochControl group="cyclone" />
              <RCPControl group="cyclone" />
            </ParamRow>
          </InputSection>
        </ToggleSection>
      </ToggleSectionGroup>
    </Box>
  );
};
