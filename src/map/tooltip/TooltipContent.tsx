import { Box, Typography } from '@mui/material';
import { FC } from 'react';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';
import { RegionHoverDescription } from './content/RegionHoverDescription';
import { useRecoilValue } from 'recoil';
import { hasHover, hoverState } from 'lib/map/interactions/interaction-state';
import { InteractionTarget } from 'lib/map/interactions/use-interactions';

export const TooltipContent: FC = () => {
  const hoveredVector = useRecoilValue(hoverState('assets')) as InteractionTarget<any>;
  const hoveredRasters = useRecoilValue(hoverState('hazards')) as InteractionTarget<any>[];
  const hoveredRegion = useRecoilValue(hoverState('regions'));
  console.log(hoveredVector, hoveredRasters, hoveredRegion);

  const assetsHovered = hasHover(hoveredVector);
  const hazardsHovered = hasHover(hoveredRasters);
  const doShow = assetsHovered || hazardsHovered;

  if (!doShow) return null;

  return (
    <>
      {assetsHovered ? (
        <Box mb={2}>
          <Typography>Asset</Typography>
          <VectorHoverDescription hoveredObject={hoveredVector} />
        </Box>
      ) : null}
      {hazardsHovered ? (
        <Box>
          <Typography>Hazards</Typography>
          {hoveredRasters.map((hr) => (
            <RasterHoverDescription hoveredObject={hr} key={`${hr.viewLayer}-${hr.target}`} />
          ))}
        </Box>
      ) : null}
      {doShow && hoveredRegion ? (
        <Box>
          <RegionHoverDescription hoveredObject={hoveredRegion} />
        </Box>
      ) : null}
    </>
  );
};
