import { Box } from '@mui/material';
import { isArray } from 'lodash';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { hasHover, hoverState } from '@/lib/data-map/interactions/interaction-state';
import { InteractionTarget } from '@/lib/data-map/interactions/use-interactions';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { AutoHidePaper, PreventHide } from './auto-hide';
import { WdpaHoverDescription } from './content/WdpaHoverDescription';

const TooltipSection = ({ children }) => (
  <Box px={1} py={0.5} borderBottom="1px solid #ccc">
    <PreventHide />
    {children}
  </Box>
);

const InteractionGroupTooltip = ({ group }) => {
  const hover = useRecoilValue(hoverState(group));

  if (!hasHover(hover)) return null;

  const contents = (
    isArray(hover)
      ? hover.map((hoveredObject) => hoveredObject.viewLayer.renderTooltip?.(hoveredObject))
      : [hover.viewLayer.renderTooltip?.(hover)]
  ).filter(Boolean);

  if (contents.length === 0) return null;

  return <TooltipSection>{contents}</TooltipSection>;
};

const WdpaTooltipSection = () => {
  const hoveredWdpas = useRecoilValue(hoverState('wdpa')) as InteractionTarget<any>[];

  if (!hasHover(hoveredWdpas)) return null;

  return (
    <TooltipSection>
      <WdpaHoverDescription hoveredObjects={hoveredWdpas} />
    </TooltipSection>
  );
};

export const TooltipContent: FC = () => {
  return (
    <AutoHidePaper>
      <Box minWidth={200}>
        <ErrorBoundary message="There was a problem displaying the tooltip.">
          <InteractionGroupTooltip group="assets" />
          <InteractionGroupTooltip group="hazards" />
          <InteractionGroupTooltip group="raster_assets" />
          <InteractionGroupTooltip group="hdi" />

          <WdpaTooltipSection />
        </ErrorBoundary>
      </Box>
    </AutoHidePaper>
  );
};
