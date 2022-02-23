import { Box, Paper, Typography } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';
import { RegionHoverDescription } from './content/RegionHoverDescription';
import { hasHover, hoverState } from 'lib/data-map/interactions/interaction-state';
import { InteractionTarget } from 'lib/data-map/interactions/use-interactions';

export const TooltipContent: FC = () => {
  const hoveredVector = useRecoilValue(hoverState('assets')) as InteractionTarget<any>;
  const hoveredRasters = useRecoilValue(hoverState('hazards')) as InteractionTarget<any>[];
  const hoveredRegion = useRecoilValue(hoverState('regions')) as InteractionTarget<any>;

  const assetsHovered = hasHover(hoveredVector);
  const hazardsHovered = hasHover(hoveredRasters);
  const doShow = assetsHovered || hazardsHovered;

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
        {doShow && hoveredRegion ? (
          <Box>
            <RegionHoverDescription hoveredObject={hoveredRegion} />
          </Box>
        ) : null}
      </Box>
    </Paper>
  );
};
