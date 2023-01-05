import { Stack } from '@mui/material';
import { Suspense } from 'react';
import { useRecoilValue } from 'recoil';

import { DataGroup } from '@/lib/data-selection/DataGroup';

import { HazardType } from '@/config/hazards/metadata';
import { DataNotice } from '@/sidebar/ui/DataNotice';
import { InputRow } from '@/sidebar/ui/InputRow';
import { EpochControl } from '@/sidebar/ui/params/EpochControl';
import { GCMControl } from '@/sidebar/ui/params/GCMControl';
import { RCPControl } from '@/sidebar/ui/params/RCPControl';
import { ReturnPeriodControl } from '@/sidebar/ui/params/ReturnPeriodControl';
import { hazardDomainsConfigState } from '@/state/data-domains/hazards';
import { paramsConfigState, useLoadParamsConfig } from '@/state/data-params';

export const InitHazardData = ({ type }: { type: HazardType }) => {
  useLoadParamsConfig(hazardDomainsConfigState(type), type);

  return null;
};

const EnsureHazardData = ({ type }) => {
  useRecoilValue(paramsConfigState(type));

  return null;
};

const HazardControl = ({ type, children }) => {
  return (
    <Suspense fallback="Loading data...">
      <EnsureHazardData type={type} />
      <DataGroup group={type}>
        <Stack spacing={3}>{children}</Stack>
      </DataGroup>
    </Suspense>
  );
};

export const FluvialControl = () => {
  return (
    <HazardControl type="fluvial">
      <DataNotice>
        Map shows river flooding depths for different return periods, from WRI Aqueduct (2020).
      </DataNotice>
      <ReturnPeriodControl />
      <InputRow>
        <EpochControl />
        <RCPControl />
      </InputRow>
    </HazardControl>
  );
};

export const CoastalControl = () => {
  return (
    <HazardControl type="coastal">
      <DataNotice>
        Map shows coastal flooding depths for different return periods, from WRI Aqueduct (2020).
      </DataNotice>
      <ReturnPeriodControl />
      <InputRow>
        <EpochControl />
        <RCPControl />
      </InputRow>
    </HazardControl>
  );
};

export const CycloneControl = () => {
  return (
    <HazardControl type="cyclone">
      <DataNotice>
        Map shows tropical cyclone maximum wind speed (in m/s) for different return periods, from
        Bloemendaal et al (2020).
      </DataNotice>
      <ReturnPeriodControl
        valueLabelDisplay="auto"
        showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
      />
      <InputRow>
        <EpochControl />
        <RCPControl />
      </InputRow>
      <GCMControl />
    </HazardControl>
  );
};

export const ExtremeHeatControl = () => {
  return (
    <HazardControl type="extreme_heat">
      <DataNotice>
        Map shows annual probability of an "extreme heat event", defined by Lange et al (2020) using
        both a relative indicator based on temperature (Russo et al 2015, 2017) and an absolute
        indicator based on temperature and relative humidity (Masterton &amp; Richardson, 1979)
        exceed their respective threshold values.
      </DataNotice>
      <InputRow>
        <EpochControl />
        <RCPControl />
      </InputRow>
      <GCMControl />
    </HazardControl>
  );
};

export const DroughtControl = () => {
  return (
    <HazardControl type="drought">
      <DataNotice>
        Map shows annual probability of a "drought event", defined by Lange et al (2020) as monthly
        soil moisture falling below the 2.5th percentile of the preindustrial baseline distribution
        for at least seven consecutive months.
      </DataNotice>
      <InputRow>
        <EpochControl />
        <RCPControl />
      </InputRow>
      <GCMControl />
    </HazardControl>
  );
};

export const EarthquakeControl = () => {
  return (
    <HazardControl type="earthquake">
      <DataNotice>
        Map shows seismic hazard as the peak ground acceleration (PGA) with a 10% probability of
        being exceeded in 50 years, from the Global Earthquake Model (GEM){' '}
        <a href="https://maps.openquake.org/map/global-seismic-hazard-map/">
          Global Seismic Hazard Map (version 2018.1)
        </a>
      </DataNotice>
    </HazardControl>
  );
};
