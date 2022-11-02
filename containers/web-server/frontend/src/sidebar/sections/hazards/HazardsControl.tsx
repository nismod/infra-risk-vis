import { Alert, Box } from '@mui/material';
import { Suspense } from 'react';
import { useRecoilValue } from 'recoil';

import { ToggleSection, ToggleSectionGroup } from '@/lib/controls/accordion-toggle/ToggleSection';

import { HazardType } from '@/config/hazards/metadata';
import { DataNotice } from '@/sidebar/ui/DataNotice';
import { InputRow } from '@/sidebar/ui/InputRow';
import { InputSection } from '@/sidebar/ui/InputSection';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { GCMControl } from '@/sidebar/ui/params/GCMControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import { ReturnPeriodControl } from '@/sidebar/ui/params/ReturnPeriodControl';
import { hazardDomainsConfigState } from '@/state/data-domains/hazards';
import { useLoadParamsConfig } from '@/state/data-params';
import { showDamagesState } from '@/state/data-selection/damage-mapping/damage-map';
import { hazardSelectionState } from '@/state/data-selection/hazards/hazard-selection';

const InitHazardData = ({ type }: { type: HazardType }) => {
  useLoadParamsConfig(hazardDomainsConfigState(type), type);

  return null;
};

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
          <Suspense fallback="Loading data...">
            <InitHazardData type="fluvial" />
            <InputSection>
              <ReturnPeriodControl group="fluvial" param="rp" disabled={disabled} />
            </InputSection>
            <InputSection>
              <InputRow>
                <EpochControl group="fluvial" disabled={disabled} />
                <RCPControl group="fluvial" disabled={disabled} />
              </InputRow>
            </InputSection>
            <InputSection>
              <GCMControl group="fluvial" disabled={disabled} />
            </InputSection>
          </Suspense>
        </ToggleSection>

        <ToggleSection id="coastal" label="Coastal Flooding" disabled={disabled}>
          <Suspense fallback="Loading data...">
            <InitHazardData type="coastal" />
            <InputSection>
              <ReturnPeriodControl group="coastal" param="rp" disabled={disabled} />
            </InputSection>
            <InputSection>
              <InputRow>
                <EpochControl group="coastal" disabled={disabled} />
                <RCPControl group="coastal" disabled={disabled} />
              </InputRow>
            </InputSection>
          </Suspense>
        </ToggleSection>

        <ToggleSection id="cyclone" label="Cyclones" disabled={disabled}>
          <Suspense fallback="Loading data...">
            <InitHazardData type="cyclone" />
            <InputSection>
              <ReturnPeriodControl
                group="cyclone"
                valueLabelDisplay="auto"
                showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
                disabled={disabled}
              />
            </InputSection>
            <InputSection>
              <GCMControl group="cyclone" disabled={disabled} />
            </InputSection>
          </Suspense>
        </ToggleSection>

        {/* <ToggleSection id="extreme_heat_exposure" label="Extreme Heat (Exposure)" disabled={disabled}>
          <InputSection>
            <InputRow>
              <EpochControl group="extreme_heat_exposure" disabled={disabled} />
              <RCPControl group="extreme_heat_exposure" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection> */}

        <ToggleSection id="extreme_heat" label="Extreme Heat" disabled={disabled}>
          <Suspense fallback="Loading data...">
            <InitHazardData type="extreme_heat" />
            <InputSection>
              <InputRow>
                <EpochControl group="extreme_heat" disabled={disabled} />
                <RCPControl group="extreme_heat" disabled={disabled} />
              </InputRow>
            </InputSection>
            <InputSection>
              <GCMControl group="extreme_heat" disabled={disabled} />
            </InputSection>
          </Suspense>
        </ToggleSection>

        <ToggleSection id="drought" label="Droughts" disabled={disabled}>
          <Suspense fallback="Loading data...">
            <InitHazardData type="drought" />
            <InputSection>
              <InputRow>
                <EpochControl group="drought" disabled={disabled} />
                <RCPControl group="drought" disabled={disabled} />
              </InputRow>
            </InputSection>
            <InputSection>
              <GCMControl group="drought" disabled={disabled} />
            </InputSection>
          </Suspense>
        </ToggleSection>

        <ToggleSection id="earthquake" label="Earthquakes">
          <Suspense fallback="Loading data...">
            <InitHazardData type="earthquake" />
            <DataNotice>
              Map shows seismic hazard as the peak ground acceleration (PGA) with a 10% probability of being exceeded in
              50 years, from the Global Earthquake Model (GEM){' '}
              <a href="https://maps.openquake.org/map/global-seismic-hazard-map/">
                Global Seismic Hazard Map (version 2018.1)
              </a>
            </DataNotice>
          </Suspense>
        </ToggleSection>

        <ToggleSection id="wildfire" label="Wildfires" disabled={true}>
          {/* Placeholder */}
        </ToggleSection>
      </ToggleSectionGroup>
    </>
  );
};
