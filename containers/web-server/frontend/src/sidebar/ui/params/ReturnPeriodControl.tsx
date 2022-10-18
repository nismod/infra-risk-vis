import { FormControl, FormLabel } from '@mui/material';

import { CustomNumberSlider } from '@/lib/controls/CustomSlider';

import { DataParam } from './DataParam';

export const ReturnPeriodControl = ({ group, ...otherProps }) => {
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
