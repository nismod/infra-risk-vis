import { Box, Paper, Typography } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';
import { RegionHoverDescription } from './content/RegionHoverDescription';
import { hasHover, hoverState } from 'lib/data-map/interactions/interaction-state';
import { InteractionTarget } from 'lib/data-map/interactions/use-interactions';
import { showPopulationState } from 'state/population';

export const TooltipContent: FC = () => {
  const hoveredVector = useRecoilValue(hoverState('assets')) as InteractionTarget<any>;
  const hoveredRasters = useRecoilValue(hoverState('hazards')) as InteractionTarget<any>[];
  const hoveredRegion = useRecoilValue(hoverState('regions')) as InteractionTarget<any>;

  const regionDataShown = useRecoilValue(showPopulationState);

  const assetsHovered = hasHover(hoveredVector);
  const hazardsHovered = hasHover(hoveredRasters);
  const regionsHovered = hasHover(hoveredRegion);
  const doShow = assetsHovered || hazardsHovered || (regionDataShown && regionsHovered);

  if (!doShow) return null;

  return (
    <Paper>
      <Box p={1}>
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
              <RasterHoverDescription hoveredObject={hr} key={`${hr.viewLayer.id}-${hr.target.id}`} />
            ))}
          </Box>
        ) : null}
        {regionsHovered ? (
          <Box>
            <RegionHoverDescription hoveredObject={hoveredRegion} />
          </Box>
        ) : null}
      </Box>
    </Paper>
  );
};
