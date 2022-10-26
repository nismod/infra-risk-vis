import { Box, Paper } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { hasHover, hoverState } from '@/lib/data-map/interactions/interaction-state';
import { InteractionTarget } from '@/lib/data-map/interactions/use-interactions';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { DroughtHoverDescription } from './content/DroughtHoverDescription';
import { HazardHoverDescription } from './content/HazardHoverDescription';
import { HdiHoverDescription } from './content/HdiHoverDescription';
import { PopulationHoverDescription } from './content/PopulationHoverDescription';
import { SolutionHoverDescription } from './content/SolutionHoverDescription';
import { VectorHoverDescription } from './content/VectorHoverDescription';

const TooltipSection = ({ children }) => (
  <Box px={1} py={0.5} borderBottom="1px solid #ccc">
    {children}
  </Box>
);

export const TooltipContent: FC = () => {
  const hoveredVector = useRecoilValue(hoverState('assets')) as InteractionTarget<any>;
  const hoveredHazards = useRecoilValue(hoverState('hazards')) as InteractionTarget<any>[];
  const hoveredPopulation = useRecoilValue(hoverState('population')) as InteractionTarget<any>;
  const hoveredRegion = useRecoilValue(hoverState('regions')) as InteractionTarget<any>;
  const hoveredSolution = useRecoilValue(hoverState('solutions')) as InteractionTarget<any>;
  const hoveredDrought = useRecoilValue(hoverState('drought')) as InteractionTarget<any>;

  const assetsHovered = hasHover(hoveredVector);
  const hazardsHovered = hasHover(hoveredHazards);
  const populationHovered = hasHover(hoveredPopulation);
  const regionsHovered = hasHover(hoveredRegion);
  const solutionsHovered = hasHover(hoveredSolution);
  const droughtHovered = hasHover(hoveredDrought);
  const doShow =
    assetsHovered || hazardsHovered || populationHovered || regionsHovered || solutionsHovered || droughtHovered;

  if (!doShow) return null;

  return (
    <Paper>
      <Box minWidth={200}>
        <ErrorBoundary message="There was a problem displaying the tooltip.">
          {solutionsHovered && (
            <TooltipSection>
              <SolutionHoverDescription hoveredObject={hoveredSolution} />
            </TooltipSection>
          )}
          {assetsHovered ? (
            <TooltipSection>
              <VectorHoverDescription hoveredObject={hoveredVector} />
            </TooltipSection>
          ) : null}
          {droughtHovered ? (
            <TooltipSection>
              <DroughtHoverDescription hoveredObject={hoveredDrought} />
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
          {regionsHovered ? (
            <TooltipSection>
              <HdiHoverDescription hoveredObject={hoveredRegion} />
            </TooltipSection>
          ) : null}
        </ErrorBoundary>
      </Box>
    </Paper>
  );
};
