import {
  FormControl,
  FormControlLabel,
  FormLabel,
  MenuItem,
  Radio,
  RadioGroup,
  Select,
} from '@mui/material';
import { useRecoilState } from 'recoil';

import { DataGroup } from '@/lib/data-selection/DataGroup';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { HAZARDS_METADATA, HAZARDS_UI_ORDER } from '@/config/hazards/metadata';
import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { LayerStylePanel } from '@/sidebar/ui/LayerStylePanel';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import {
  damageSourceState,
  damageSourceStateEffect,
} from '@/state/data-selection/damage-mapping/damage-map';

export const DamageSourceControl = () => {
  const [damageSource, setDamageSource] = useRecoilState(damageSourceState);

  return (
    <>
      <StateEffectRoot state={damageSourceState} effect={damageSourceStateEffect} />
      <LayerStylePanel>
        <InputSection>
          <FormControl>
            <FormLabel>Hazard</FormLabel>
            <RadioGroup value={damageSource} onChange={(e, value) => setDamageSource(value)}>
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
          <DataGroup group={damageSource}>
            <InputRow>
              <EpochControl />
              <RCPControl />
            </InputRow>
          </DataGroup>
        </InputSection>
      </LayerStylePanel>
    </>
  );
};
