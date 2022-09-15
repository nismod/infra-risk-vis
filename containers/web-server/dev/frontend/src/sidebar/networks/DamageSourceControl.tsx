import {
  Box,
  FormControl,
  FormControlLabel,
  FormLabel,
  MenuItem,
  Paper,
  Radio,
  RadioGroup,
  Select,
} from '@mui/material';
import { useRecoilState } from 'recoil';

import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';

import { InputSection } from 'sidebar/ui/InputSection';
import { InputRow } from 'sidebar/ui/InputRow';
import { EpochControl } from 'sidebar/ui/params/EpochControl';
import { RCPControl } from 'sidebar/ui/params/RCPControl';
import { damageSourceState, damageSourceStateEffect, damageTypeState } from 'state/damage-mapping/damage-map';
import { HAZARDS_METADATA, HAZARDS_UI_ORDER } from 'config/hazards/metadata';
import { LayerStylePanel } from 'sidebar/ui/LayerStylePanel';

export const DamageSourceControl = () => {
  const [damageSource, setDamageSource] = useRecoilState(damageSourceState);
  const [damageType, setDamageType] = useRecoilState(damageTypeState);

  return (
    <>
      <StateEffectRoot state={damageSourceState} effect={damageSourceStateEffect} />
      <LayerStylePanel>
        <InputSection>
          <FormControl fullWidth>
            <FormLabel>Damage type</FormLabel>
            <Select<string> variant="standard" value={damageType} onChange={(e) => setDamageType(e.target.value)}>
              <MenuItem value="direct">Direct Damages</MenuItem>
              <MenuItem value="indirect">Economic Losses</MenuItem>
            </Select>
          </FormControl>
        </InputSection>
        <InputSection>
          <FormControl>
            <FormLabel>Hazard</FormLabel>
            <RadioGroup value={damageSource} onChange={(e, value) => setDamageSource(value)}>
              <FormControlLabel label="All Hazards" control={<Radio value="all" />} />
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
      </LayerStylePanel>
    </>
  );
};
