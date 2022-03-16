import { Box, FormControl, FormControlLabel, FormLabel, Paper, Radio, RadioGroup } from '@mui/material';
import { useRecoilState } from 'recoil';

import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';

import { InputSection } from 'sidebar/ui/InputSection';
import { InputRow } from 'sidebar/ui/InputRow';
import { EpochControl } from 'sidebar/ui/params/EpochControl';
import { RCPControl } from 'sidebar/ui/params/RCPControl';
import { damageSourceState, damageSourceStateEffect } from 'state/damage-mapping/damage-map';
import { HAZARDS_METADATA, HAZARDS_UI_ORDER } from 'config/hazards/metadata';

export const DamageSourceControl = () => {
  const [damageSource, setDamageSource] = useRecoilState(damageSourceState);

  return (
    <>
      <StateEffectRoot state={damageSourceState} effect={damageSourceStateEffect} />
      <Box mt={2}>
        <Paper>
          <Box p={2}>
            <InputSection>
              <FormControl>
                <FormLabel>Hazard</FormLabel>
                <RadioGroup value={damageSource} onChange={(e, value) => setDamageSource(value)}>
                  <FormControlLabel label="Total Damages" control={<Radio value="total-damages" />} />
                  {HAZARDS_UI_ORDER.map((hazard) => (
                    <FormControlLabel
                      key={hazard}
                      label={HAZARDS_METADATA[hazard].label}
                      control={<Radio value={hazard} />}
                    />
                  ))}
                </RadioGroup>
              </FormControl>
            </InputSection>
            <InputSection>
              <InputRow>
                <EpochControl group={damageSource} />
                <RCPControl group={damageSource} />
              </InputRow>
            </InputSection>
          </Box>
        </Paper>
      </Box>
    </>
  );
};
