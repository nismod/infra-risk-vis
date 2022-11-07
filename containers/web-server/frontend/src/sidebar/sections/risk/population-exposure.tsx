import { FormControl, FormControlLabel, FormLabel, Radio, RadioGroup, Switch } from '@mui/material';
import _ from 'lodash';
import { useEffect } from 'react';
import { TransactionInterface_UNSTABLE, atom, useRecoilState, useRecoilTransaction_UNSTABLE } from 'recoil';

import { DataGroup } from '@/lib/data-selection/DataGroup';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { ExposureSource } from '@/config/hazards/exposure/exposure-view-layer';
import { sidebarPathChildrenState, sidebarVisibilityToggleState } from '@/sidebar/SidebarContent';
import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import { syncHazardsWithDamageSourceStateEffect } from '@/state/data-selection/damage-mapping/damage-map';

export const populationExposureHazardState = atom<ExposureSource>({
  key: 'populationExposureHazardState',
  default: 'extreme_heat',
});

/**
 * set only population layer to visible
 */
function syncExposure({ get, set }: TransactionInterface_UNSTABLE) {
  const hazardSubPaths = get(sidebarPathChildrenState('exposure'));

  _.forEach(hazardSubPaths, (subPath) => {
    set(sidebarVisibilityToggleState(`exposure/${subPath}`), subPath === 'population');
  });
}

export const PopulationExposureSection = () => {
  const [hazard, setHazard] = useRecoilState(populationExposureHazardState);

  const [showHazards, setShowHazards] = useRecoilState(sidebarVisibilityToggleState('hazards'));
  const [showPopulation, setShowPopulation] = useRecoilState(sidebarVisibilityToggleState('exposure'));

  const updateExposure = useRecoilTransaction_UNSTABLE((iface) => () => syncExposure(iface));

  useEffect(() => {
    updateExposure();
  }, [updateExposure]);

  return (
    <>
      <StateEffectRoot state={populationExposureHazardState} effect={syncHazardsWithDamageSourceStateEffect} />
      <InputSection>
        <FormControl>
          <FormLabel>Hazard</FormLabel>
          <RadioGroup value={hazard} onChange={(e, value) => setHazard(value as ExposureSource)}>
            <FormControlLabel value="extreme_heat" control={<Radio />} label="Extreme Heat" />
            <FormControlLabel value="drought" control={<Radio />} label="Droughts" />
          </RadioGroup>
        </FormControl>
      </InputSection>
      <InputSection>
        <DataGroup group={hazard}>
          <InputRow>
            <EpochControl />
            <RCPControl />
          </InputRow>
        </DataGroup>
      </InputSection>
      <InputSection>
        <FormControlLabel
          control={<Switch />}
          value={showHazards}
          onChange={(e, checked) => setShowHazards(checked)}
          label="Hazard layer"
        />
      </InputSection>
      <InputSection>
        <FormControlLabel
          control={<Switch />}
          value={showPopulation}
          onChange={(e, checked) => setShowPopulation(checked)}
          label="Population layer"
        />
      </InputSection>
    </>
  );
};
