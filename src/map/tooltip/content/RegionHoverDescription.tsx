import { useRecoilValue } from 'recoil';
import { Typography } from '@mui/material';

import { boundaryConfig, BoundaryLevel } from '../../../config/deck-layers/boundaries-layer';
import { boundaryLevelState } from '../../layers/layers-state';
import { InteractionTarget, VectorTarget } from 'lib/map/interactions/use-interactions';

const labels: Record<BoundaryLevel, string> = {
  parish: 'Parish',
  community: 'Community',
  enumeration: 'Enumeration District',
};

export const RegionHoverDescription = ({ hoveredObject }: { hoveredObject: InteractionTarget<VectorTarget> }) => {
  const boundaryLevel = useRecoilValue(boundaryLevelState);
  const levelLabel = labels[boundaryLevel];

  return (
    <>
      <Typography component="h6">{levelLabel}</Typography>
      <Typography>{hoveredObject.target.feature.properties[boundaryConfig[boundaryLevel].fieldName]}</Typography>
    </>
  );
};
