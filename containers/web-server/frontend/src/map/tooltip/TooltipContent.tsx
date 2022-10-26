import { Box, Paper } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { hasHover, hoverState } from '@/lib/data-map/interactions/interaction-state';
import { InteractionTarget } from '@/lib/data-map/interactions/use-interactions';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { NETWORKS_METADATA } from '@/config/networks/metadata';

import { VectorHoverDescription } from './VectorHoverDescription';
import { HazardHoverDescription } from './content/HazardHoverDescription';
import { HdiHoverDescription } from './content/HdiHoverDescription';
import { PopulationHoverDescription } from './content/PopulationHoverDescription';
import { WdpaHoverDescription } from './content/WdpaHoverDescription';

const TooltipSection = ({ children }) => (
  <Box px={1} py={0.5} borderBottom="1px solid #ccc">
    {children}
  </Box>
);

export const TooltipContent: FC = () => {
  const hoveredVector = useRecoilValue(hoverState('assets')) as InteractionTarget<any>;
  const hoveredHazards = useRecoilValue(hoverState('hazards')) as InteractionTarget<any>[];
  const hoveredPopulation = useRecoilValue(hoverState('population')) as InteractionTarget<any>;
  const hoveredHdi = useRecoilValue(hoverState('hdi')) as InteractionTarget<any>;
  const hoveredWdpas = useRecoilValue(hoverState('wdpa')) as InteractionTarget<any>[];

  const assetsHovered = hasHover(hoveredVector);
  const hazardsHovered = hasHover(hoveredHazards);
  const populationHovered = hasHover(hoveredPopulation);
  const hdiHovered = hasHover(hoveredHdi);
  const wdpaHovered = hasHover(hoveredWdpas);

  const doShow = assetsHovered || hazardsHovered || populationHovered || hdiHovered || wdpaHovered;

  if (!doShow) return null;

  return (
    <Paper>
      <Box minWidth={200}>
        <ErrorBoundary message="There was a problem displaying the tooltip.">
          {/* TODO: generate tooltip contents straight from view layers */}
          {assetsHovered ? (
            <TooltipSection>
              <VectorHoverDescription hoveredObject={hoveredVector} metadataLookup={NETWORKS_METADATA} />
            </TooltipSection>
          ) : null}
          {hazardsHovered ? (
            <TooltipSection>
              {hoveredHazards.map((hr) => (
                <HazardHoverDescription hoveredObject={hr} key={`${hr.viewLayer.id}-${hr.target.id}`} />
              ))}
            </TooltipSection>
          ) : null}
          {populationHovered ? (
            <TooltipSection>
              <PopulationHoverDescription hoveredObject={hoveredPopulation} />
            </TooltipSection>
          ) : null}
          {hdiHovered ? (
            <TooltipSection>
              <HdiHoverDescription hoveredObject={hoveredHdi} />
            </TooltipSection>
          ) : null}
          {wdpaHovered ? (
            <TooltipSection>
              <WdpaHoverDescription hoveredObjects={hoveredWdpas} />
            </TooltipSection>
          ) : null}
        </ErrorBoundary>
      </Box>
    </Paper>
  );
};
