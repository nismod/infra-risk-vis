import { REGIONS_METADATA } from '../../../config/regions/metadata';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { useRecoilValue } from 'recoil';
import { showPopulationState } from 'state/regions';
import { DataItem } from 'details/features/detail-components';

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
