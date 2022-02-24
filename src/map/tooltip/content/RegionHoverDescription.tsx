import { Typography } from '@mui/material';

import { REGIONS_METADATA } from '../../../config/regions/metadata';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { useRecoilValue } from 'recoil';
import { showPopulationState } from 'state/regions';

export const RegionHoverDescription = ({ hoveredObject }: { hoveredObject: InteractionTarget<VectorTarget> }) => {
  const metadata = REGIONS_METADATA[hoveredObject.viewLayer.params.boundaryLevel];

  const showPopulation = useRecoilValue(showPopulationState);

  return (
    <>
      <Typography component="h6">{metadata.labelSingular}</Typography>
      <Typography>{hoveredObject.target.feature.properties[metadata.fieldName]}</Typography>
      {showPopulation && (
        <Typography>Population: {hoveredObject.target.feature.properties['population'].toLocaleString()}</Typography>
      )}
    </>
  );
};
