import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';
import { useCallback } from 'react';

export const ParamDropdown = ({ title, value, onChange, options, variant = undefined }) => {
  const handleChange = useCallback((e) => onChange(e.target.value), [onChange]);
  return (
    <FormControl fullWidth>
      <FormLabel>{title}</FormLabel>
      <Select value={value} onChange={handleChange} size="small" variant={variant} disabled={options.length < 2}>
        {options.map((option) => {
          let value, label;
          if (typeof option === 'string' || typeof option === 'number') {
            value = label = option;
          } else {
            value = option.value;
            label = option.label;
          }

          return (
            <MenuItem key={value} value={value}>
              {label}
            </MenuItem>
          );
        })}
      </Select>
    </FormControl>
  );
};
