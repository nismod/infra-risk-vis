import { FormControl, FormLabel } from '@mui/material';

import { CustomNumberSlider } from '@/lib/controls/CustomSlider';

import { useDataGroup } from '../../../lib/data-selection/DataGroup';
import { DataParam } from '../../../lib/data-selection/DataParam';

export const ReturnPeriodControl = ({ ...otherProps }) => {
  const group = useDataGroup();
  return (
    <FormControl fullWidth>
      <FormLabel>Return Period</FormLabel>
      <DataParam group={group} id="rp">
        {({ value, onChange, options }) => (
          <CustomNumberSlider marks={options} value={value} onChange={onChange} size="small" {...otherProps} />
        )}
      </DataParam>
    </FormControl>
  );
};
