import { Box, Typography } from '@mui/material';
import { FC } from 'react';

import { InteractionTarget } from '@/lib/data-map/interactions/use-interactions';
import { ColorBox } from '@/lib/ui/data-display/ColorBox';

import { PROTECTED_AREA_COLORS } from '@/config/protected-areas/metadata';

export const WdpaHoverDescription: FC<{
  hoveredObjects: InteractionTarget<any>[];
}> = ({ hoveredObjects }) => {
  return (
    <>
      <Typography variant="body2">Protected Areas</Typography>
      {hoveredObjects.map((ho) => {
        const feature = ho.target.feature;
        const color = PROTECTED_AREA_COLORS[ho.viewLayer.params.type].css;
        return (
          <Box>
            <ColorBox color={color} />
            {feature.properties.NAME}
          </Box>
        );
      })}
    </>
  );
};
