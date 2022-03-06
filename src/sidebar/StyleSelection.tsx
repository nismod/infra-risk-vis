import { FormControl, InputLabel, MenuItem, Select } from '@mui/material';
import { uniqueId } from 'lodash';
import { FC, useState } from 'react';
import { RecoilState, useRecoilState } from 'recoil';

interface StyleSelectionOption {
  label: string;
  value: string;
}

export const StyleSelection: FC<{ state: RecoilState<string>; options: StyleSelectionOption[]; defaultValue: string }> =
  ({ state, options }) => {
    const id = useState(uniqueId('style-selection-'));
    const labelId = `${id}-input-label`;

    const [value, setValue] = useRecoilState(state);

    return (
      <FormControl fullWidth>
        <InputLabel id={labelId}>Layer style</InputLabel>
        <Select labelId={labelId} label="Layer style" value={value} onChange={(e) => setValue(e.target.value)}>
          {options.map((option) => (
            <MenuItem key={option.value} value={option.value}>
              {option.label}
            </MenuItem>
          ))}
        </Select>
      </FormControl>
    );
  };
