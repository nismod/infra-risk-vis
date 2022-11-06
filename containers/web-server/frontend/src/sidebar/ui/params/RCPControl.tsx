import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';

import { useDataGroup } from '@/lib/data-selection/DataGroup';
import { DataParam } from '@/lib/data-selection/DataParam';
import { useInputDisabled } from '@/lib/data-selection/DisabledInput';

function rcpLabel(value) {
  return value === 'baseline' ? 'Baseline' : value;
}

export const RCPControl = () => {
  const group = useDataGroup();
  const disabled = useInputDisabled();

  return (
    <FormControl fullWidth disabled={disabled}>
      <FormLabel>
        <abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr>
      </FormLabel>
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
