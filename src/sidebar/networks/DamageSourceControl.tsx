import { Box, FormControl, FormControlLabel, FormLabel, Paper, Radio, RadioGroup } from '@mui/material';
import { useRecoilValue } from 'recoil';
import { InputSection } from 'sidebar/ui/InputSection';
import { InputRow } from 'sidebar/ui/InputRow';
import { EpochControl } from 'sidebar/ui/params/EpochControl';
import { RCPControl } from 'sidebar/ui/params/RCPControl';

import { selectedDamageSourceState, useUpdateDamageSource } from 'state/damage-mapping/damage-map';
import { useInitDirectDamages } from '../../state/damage-mapping/damage-map';
import { useEffect } from 'react';

export const DamageSourceControl = () => {
  const selectedDamageSource = useRecoilValue(selectedDamageSourceState);
  const updateDamageSource = useUpdateDamageSource();

  const initDirectDamages = useInitDirectDamages();

  useEffect(() => {
    initDirectDamages();
  }, [initDirectDamages]);

  return (
    <Box my={2}>
      <Paper>
        <Box p={2}>
          <InputSection>
            <FormControl>
              <FormLabel>Damage Source</FormLabel>
              <RadioGroup value={selectedDamageSource} onChange={(e, value) => updateDamageSource(value)}>
                <FormControlLabel label="Total Damages" control={<Radio value="total-damages" />} />
                <FormControlLabel label="River Flooding" control={<Radio value="fluvial" />} />
                <FormControlLabel label="Surface Flooding" control={<Radio value="surface" />} />
                <FormControlLabel label="Coastal Flooding" control={<Radio value="coastal" />} />
                <FormControlLabel label="Cyclones" control={<Radio value="cyclone" />} />
              </RadioGroup>
            </FormControl>
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group={selectedDamageSource} />
              <RCPControl group={selectedDamageSource} />
            </InputRow>
          </InputSection>
        </Box>
      </Paper>
    </Box>
  );
};
