import { useRecoilValue } from 'recoil';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { REGIONS_METADATA } from '@/config/_old/regions/metadata';
import { showPopulationState } from '@/state/data-selection/_old/regions';

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
