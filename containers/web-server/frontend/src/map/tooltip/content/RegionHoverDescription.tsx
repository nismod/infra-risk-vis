import { useRecoilValue } from 'recoil';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';

import { REGIONS_METADATA } from '@/config/regions/metadata';
import { DataItem } from '@/details/features/detail-components';
import { showPopulationState } from '@/state/regions';

export const RegionHoverDescription = ({ hoveredObject }: { hoveredObject: InteractionTarget<VectorTarget> }) => {
  const metadata = REGIONS_METADATA[hoveredObject.viewLayer.params.regionLevel];

  const showPopulation = useRecoilValue(showPopulationState);

  return (
    <>
      <DataItem label={metadata.labelSingular} value={hoveredObject.target.feature.properties[metadata.fieldName]} />
      {showPopulation && (
        <DataItem label="Population" value={hoveredObject.target.feature.properties.population.toLocaleString()} />
      )}
    </>
  );
};
