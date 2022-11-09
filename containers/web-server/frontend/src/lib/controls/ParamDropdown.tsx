import { FormControl, FormLabel, MenuItem, Select, SelectProps } from '@mui/material';
import { PropsWithChildren, useCallback } from 'react';

import { ValueLabel, isValueLabel } from './params/value-label';

interface ParamDropdownProps<V extends string | number = string> {
  title?: string;
  value: V;
  options: (V | ValueLabel<V>)[];
  onChange: (value: V) => void;
  disabled?: boolean;
  variant?: SelectProps['variant'];
}

export const ParamDropdown = <V extends string | number = string>({
  title,
  value,
  onChange,
  options,
  disabled = false,
  variant = undefined,
}: PropsWithChildren<ParamDropdownProps<V>>) => {
  const handleChange = useCallback((e) => onChange(e.target.value), [onChange]);
  return (
    <FormControl fullWidth>
      {title && <FormLabel>{title}</FormLabel>}
      <Select
        value={value}
        onChange={handleChange}
        size="small"
        variant={variant}
        disabled={disabled || options.length < 2}
      >
        {options.map((option) => {
          let value, label;
          if (isValueLabel(option)) {
            ({ value, label } = option);
          } else {
            value = label = option;
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
