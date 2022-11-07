import { FC } from 'react';

import { Layer } from '@/lib/data-selection/sidebar/components';

import { PopulationExposureSection } from './population-exposure';

export const RiskSection: FC = () => {
  return (
    <>
      <Layer path="population" title="Population Exposure">
        <PopulationExposureSection />
      </Layer>
      <Layer path="infrastructure" title="Infrastructure Risk"></Layer>
      <Layer path="regional" title="Regional Risk"></Layer>
    </>
  );
};
