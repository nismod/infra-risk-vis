import { Box, Typography } from '@mui/material';
import { FC } from 'react';

import { RasterHover, RegionHover, VectorHover } from '../DataMap';

import { VectorHoverDescription } from './content/VectorHoverDescription';
import { RasterHoverDescription } from './content/RasterHoverDescription';
import { RegionHoverDescription } from './content/RegionHoverDescription';

export const TooltipContent: FC<{
  hoveredVectors: VectorHover[];
  hoveredRasters: RasterHover[];
  hoveredRegion: RegionHover;
}> = ({ hoveredVectors, hoveredRasters, hoveredRegion }) => {
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
      {(hoveredRasters.length || hoveredVectors.length) && hoveredRegion ? (
        <Box>
          <RegionHoverDescription hoveredObject={hoveredRegion} />
        </Box>
      ) : null}
    </>
  );
};
