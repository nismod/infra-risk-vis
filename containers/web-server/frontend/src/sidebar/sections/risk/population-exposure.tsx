import { FormControl, FormControlLabel, FormLabel, Radio, RadioGroup, Switch } from '@mui/material';
import _ from 'lodash';
import { useEffect } from 'react';
import { TransactionInterface_UNSTABLE, atom, useRecoilState, useRecoilTransaction_UNSTABLE } from 'recoil';

import { DataGroup } from '@/lib/data-selection/DataGroup';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { ExposureSource } from '@/config/hazards/exposure/exposure-view-layer';
import {
  sidebarPathChildrenState,
  sidebarPathVisibilityState,
  sidebarVisibilityToggleState,
} from '@/sidebar/SidebarContent';
import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import { syncHazardsWithDamageSourceStateEffect } from '@/state/data-selection/damage-mapping/damage-map';
import { DataNotice } from '@/sidebar/ui/DataNotice';

export const populationExposureHazardState = atom<ExposureSource>({
  key: 'populationExposureHazardState',
  default: 'extreme_heat',
});

/**
 * set only population layer to visible
 */
export function syncExposure({ get, set }: TransactionInterface_UNSTABLE, layer: string) {
  const hazardSubPaths = get(sidebarPathChildrenState('exposure'));

  /**
   * Using sidebarVisibilityToggleState here for individual levels
   * instead of the recursive sidebarPathVisibilityState
   * because Recoil doesn't allow setting selectors in atomic transactions
   */

  set(sidebarVisibilityToggleState('exposure'), true);
  _.forEach(hazardSubPaths, (subPath) => {
    set(sidebarVisibilityToggleState(`exposure/${subPath}`), subPath === layer);
  });
}

export function hideExposure({ set }: TransactionInterface_UNSTABLE, layer: string) {
  set(sidebarVisibilityToggleState(`exposure/${layer}`), false);
}

const InitPopulationView = () => {
  const updateExposureTx = useRecoilTransaction_UNSTABLE((iface) => () => syncExposure(iface, 'population'), []);
  const hideExposureTx = useRecoilTransaction_UNSTABLE((iface) => () => hideExposure(iface, 'population'), []);
  useEffect(() => {
    updateExposureTx();

    return () => {
      hideExposureTx();
    };
  }, [updateExposureTx, hideExposureTx]);

  return null;
};

export const PopulationExposureSection = () => {
  const [hazard, setHazard] = useRecoilState(populationExposureHazardState);

  const [showHazards, setShowHazards] = useRecoilState(sidebarPathVisibilityState(`hazards/${hazard}`));
  const [showPopulation, setShowPopulation] = useRecoilState(sidebarPathVisibilityState('exposure/population'));

  return (
    <>
      <InitPopulationView />
      <StateEffectRoot state={populationExposureHazardState} effect={syncHazardsWithDamageSourceStateEffect} />
      <DataNotice>
        Map shows expected annual population exposed to extreme events, based
        on the annual probability of the hazard. Zoom in for extreme heat
        exposure to show.<br/>
      </DataNotice>
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
          checked={showHazards}
          onChange={(e, checked) => setShowHazards(checked)}
          label="Hazard layer"
        />
      </InputSection>
      <InputSection>
        <FormControlLabel
          control={<Switch />}
          checked={showPopulation}
          onChange={(e, checked) => setShowPopulation(checked)}
          label="Population layer"
        />
      </InputSection>
    </>
  );
};
