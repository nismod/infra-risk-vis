import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';

import { DataParam } from './DataParam';

function epochLabel(value) {
  if (value === 2010) return 'Present';
  return value;
}

export const EpochControl = ({ group, disabled = false }) => {
  return (
    <FormControl fullWidth disabled={disabled}>
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
