import { Box, Typography } from '@material-ui/core';
import { FC } from 'react';

import { RasterHover, VectorHover } from '../DataMap';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';

export const TooltipContent: FC<{
  hoveredVectors: VectorHover[];
  hoveredRasters: RasterHover[];
}> = ({ hoveredVectors, hoveredRasters }) => {
  return (
    <>
      {hoveredVectors.length ? (
        <Box mb={2}>
          <Typography>Asset</Typography>
          {hoveredVectors.map((hv) => (
            <VectorHoverDescription hoveredObject={hv} key={hv.feature.id} />
          ))}
        </Box>
      ) : null}
      {hoveredRasters.length ? (
        <Box>
          <Typography>Hazards</Typography>
          {hoveredRasters.map((hr) => (
            <RasterHoverDescription hoveredObject={hr} key={`${hr.deckLayer}-${hr.color}`} />
          ))}
        </Box>
      ) : null}
    </>
  );
};
