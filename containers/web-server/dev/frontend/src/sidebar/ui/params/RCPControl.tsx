import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';
import { DataParam } from './DataParam';

function rcpLabel(value) {
  return value === 'baseline' ? 'Baseline' : value;
}

export const RCPControl = ({ group, disabled = false }) => {
  return (
    <FormControl fullWidth disabled={disabled}>
      <FormLabel><abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr></FormLabel>
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
