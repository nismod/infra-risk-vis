import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';

import { useDataGroup } from '../../../lib/data-selection/DataGroup';
import { DataParam } from '../../../lib/data-selection/DataParam';
import { useInputDisabled } from '../../../lib/data-selection/DisabledInput';

function epochLabel(value) {
  if (value === 2010 || value === 'present' || value === 'baseline') return 'Present';
  return value;
}

export const EpochControl = () => {
  const group = useDataGroup();
  const disabled = useInputDisabled();
  return (
    <FormControl fullWidth disabled={disabled}>
      <FormLabel>Epoch</FormLabel>
      <DataParam group={group} id="epoch">
        {({ value, onChange, options }) => (
          <Select
            variant="standard"
            value={value}
            onChange={(e) => onChange(e.target.value)}
            fullWidth
          >
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
