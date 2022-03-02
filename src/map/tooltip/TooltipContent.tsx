import { Box, Paper } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';
import { RegionHoverDescription } from './content/RegionHoverDescription';
import { hasHover, hoverState } from 'lib/data-map/interactions/interaction-state';
import { InteractionTarget } from 'lib/data-map/interactions/use-interactions';
import { showPopulationState } from 'state/regions';

const TooltipSection = ({ children }) => (
  <Box p={1} borderBottom="1px solid #ccc">
    {children}
  </Box>
);

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
      <Box minWidth={200}>
        {assetsHovered ? (
          <TooltipSection>
            <VectorHoverDescription hoveredObject={hoveredVector} />
          </TooltipSection>
        ) : null}
        {hazardsHovered ? (
          <TooltipSection>
            {hoveredRasters.map((hr) => (
              <RasterHoverDescription hoveredObject={hr} key={`${hr.viewLayer.id}-${hr.target.id}`} />
            ))}
          </TooltipSection>
        ) : null}
        {regionsHovered ? (
          <TooltipSection>
            <RegionHoverDescription hoveredObject={hoveredRegion} />
          </TooltipSection>
        ) : null}
      </Box>
    </Paper>
  );
};
