import { Typography } from '@mui/material';

import { REGIONS_METADATA } from '../../../config/regions/metadata';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';

export const RegionHoverDescription = ({ hoveredObject }: { hoveredObject: InteractionTarget<VectorTarget> }) => {
  const metadata = REGIONS_METADATA[hoveredObject.viewLayer.params.boundaryLevel];

  return (
    <>
      <Typography component="h6">{metadata.label}</Typography>
      <Typography>{hoveredObject.target.feature.properties[metadata.fieldName]}</Typography>
    </>
  );
};
