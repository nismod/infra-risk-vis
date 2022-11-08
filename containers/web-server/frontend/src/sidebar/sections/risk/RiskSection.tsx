import { FC } from 'react';

import { Layer } from '@/lib/data-selection/sidebar/components';
import { EnforceSingleChild } from '@/lib/data-selection/sidebar/single-child';

import { InfrastructureRiskSection } from './infrastructure-risk';
import { PopulationExposureSection } from './population-exposure';

export const RiskSection: FC = () => {
  return (
    <>
      <EnforceSingleChild />
      <Layer path="population" title="Population Exposure" unmountOnHide={true}>
        <PopulationExposureSection />
      </Layer>
      <Layer path="infrastructure" title="Infrastructure Risk" unmountOnHide={true}>
        <InfrastructureRiskSection />
      </Layer>
      <Layer path="regional" title="Regional Risk" disabled></Layer>
    </>
  );
};
