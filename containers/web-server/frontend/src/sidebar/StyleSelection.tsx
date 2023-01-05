import { FormControl, InputLabel, MenuItem, Select } from '@mui/material';
import { uniqueId } from 'lodash';
import { FC, useState } from 'react';
import { useRecoilState, useRecoilValue } from 'recoil';

import { sectionStyleOptionsState, sectionStyleValueState } from '@/state/sections';

export const StyleSelection: FC<{ id: string }> = ({ id }) => {
  const [value, setValue] = useRecoilState(sectionStyleValueState(id));
  const options = useRecoilValue(sectionStyleOptionsState(id));

  const htmlId = useState(uniqueId('style-selection-'));
  const labelId = `${htmlId}-input-label`;

  return (
    <FormControl fullWidth>
      <InputLabel id={labelId}>Layer style</InputLabel>
      <Select
        labelId={labelId}
        label="Layer style"
        value={value}
        onChange={(e) => setValue(e.target.value)}
      >
        {options.map((option) => (
          <MenuItem key={option.id} value={option.id}>
            {option.label}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  );
};
