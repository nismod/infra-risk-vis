import { FormControl, FormLabel, MenuItem, Select } from '@mui/material';

import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { DataParam } from '@/sidebar/ui/params/DataParam';

export const EarthquakeToggleSection = ({ disabled = false }) => {
  return (
    <InputSection>
      <InputRow>
        <FormControl fullWidth disabled={disabled}>
          <FormLabel>Return Period</FormLabel>
          <DataParam group="earthquake" id="returnPeriod">
            {({ value, onChange, options }) => {
              return (
                <Select variant="standard" value={value} onChange={(e) => onChange(e.target.value)} fullWidth>
                  {options.map((rp) => (
                    <MenuItem key={rp} value={rp}>
                      {rp}
                    </MenuItem>
                  ))}
                </Select>
              );
            }}
          </DataParam>
        </FormControl>
        <FormControl fullWidth disabled={disabled}>
          <FormLabel>Medium</FormLabel>
          <DataParam group="earthquake" id="medium">
            {({ value, onChange, options }) => {
              return (
                <Select variant="standard" value={value} onChange={(e) => onChange(e.target.value)} fullWidth>
                  {options.map((rp) => (
                    <MenuItem key={rp} value={rp}>
                      {rp}
                    </MenuItem>
                  ))}
                </Select>
              );
            }}
          </DataParam>
        </FormControl>
      </InputRow>
    </InputSection>
  );
};
