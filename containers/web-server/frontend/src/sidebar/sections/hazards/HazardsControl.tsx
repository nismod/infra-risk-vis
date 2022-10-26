import { Alert, Box } from '@mui/material';
import { useRecoilValue } from 'recoil';

import { ToggleSection, ToggleSectionGroup } from '@/lib/controls/accordion-toggle/ToggleSection';

import { DataNotice } from '@/sidebar/ui/DataNotice';
import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import { ReturnPeriodControl } from '@/sidebar/ui/params/ReturnPeriodControl';
import { showDamagesState } from '@/state/damage-mapping/damage-map';
import { hazardSelectionState } from '@/state/hazards/hazard-selection';

export const HazardsControl = () => {
  const showDirectDamages = useRecoilValue(showDamagesState);
  const disabled = showDirectDamages;

  return (
    <>
      {showDirectDamages ? (
        <Box my={1}>
          <Alert severity="info">Hazards are currently following the Infrastructure &gt; Damages selection</Alert>
        </Box>
      ) : null}
      <ToggleSectionGroup toggleState={hazardSelectionState}>
        <ToggleSection id="fluvial" label="River Flooding" disabled={disabled}>
          <InputSection>
            <ReturnPeriodControl group="fluvial" param="returnPeriod" disabled={disabled} />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="fluvial" disabled={disabled} />
              <RCPControl group="fluvial" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="coastal" label="Coastal Flooding" disabled={disabled}>
          <InputSection>
            <ReturnPeriodControl group="coastal" param="returnPeriod" disabled={disabled} />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="coastal" disabled={disabled} />
              <RCPControl group="coastal" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="cyclone" label="Cyclones" disabled={disabled}>
          <InputSection>
            <ReturnPeriodControl
              group="cyclone"
              valueLabelDisplay="auto"
              showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
              disabled={disabled}
            />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="cyclone" disabled={disabled} />
              <RCPControl group="cyclone" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="extreme_heat_exposure" label="Extreme Heat (Exposure)" disabled={disabled}>
          <InputSection>
            <InputRow>
              <EpochControl group="extreme_heat_exposure" disabled={disabled} />
              <RCPControl group="extreme_heat_exposure" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="extreme_heat_occurrence" label="Extreme Heat (Occurrence)" disabled={disabled}>
          <InputSection>
            <InputRow>
              <EpochControl group="extreme_heat_occurrence" disabled={disabled} />
              <RCPControl group="extreme_heat_occurrence" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>
        <ToggleSection id="drought" label="Droughts" disabled={true}>
          {/* Placeholder */}
        </ToggleSection>
        <ToggleSection id="wildfire" label="Wildfires" disabled={true}>
          {/* Placeholder */}
        </ToggleSection>
        <ToggleSection id="earthquake" label="Earthquakes">
          <DataNotice>
            Map shows seismic hazard as the peak ground acceleration (PGA) with a 10% probability of being exceeded in
            50 years, from the Global Earthquake Model (GEM){' '}
            <a href="https://maps.openquake.org/map/global-seismic-hazard-map/">
              Global Seismic Hazard Map (version 2018.1)
            </a>
          </DataNotice>
        </ToggleSection>
      </ToggleSectionGroup>
    </>
  );
};
