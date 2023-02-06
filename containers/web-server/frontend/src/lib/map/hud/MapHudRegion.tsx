import { PropsOf } from '@emotion/react';
import { Box, Stack } from '@mui/material';
import React, { FC } from 'react';

import { eventPreventDefault } from '@/lib/helpers';

const hudRegions = {
  'top-left': {
    style: {
      top: 0,
      left: 0,
    },
    justifyContent: 'left',
  },
  'top-right': {
    style: {
      top: 0,
      right: 0,
    },
    justifyContent: 'right',
  },
  'bottom-right': {
    style: {
      bottom: 0,
      right: 0,
    },
    justifyContent: 'right',
  },
  'bottom-left': {
    style: {
      bottom: 0,
      left: 0,
    },
    justifyContent: 'left',
  },
};

export interface MapHudProps {
  position: keyof typeof hudRegions;
  style?: PropsOf<typeof Box>['style'];
  StackProps?: PropsOf<typeof Stack>;
}
export const MapHudRegion: FC<MapHudProps> = ({ position, style = {}, StackProps, children }) => {
  const { style: defaultRegionStyle, justifyContent } = hudRegions[position];
  const effectiveStyle = {
    ...defaultRegionStyle,
    m: 1,
    style,
  };
  return (
    <Box
      position="absolute"
      {...effectiveStyle}
      // stop interactions with HUD from triggering map events: https://github.com/visgl/deck.gl/discussions/6252
      onPointerUp={eventPreventDefault}
      onMouseMove={eventPreventDefault}
      sx={{ pointerEvents: 'auto' }}
    >
      <Stack {...StackProps}>
        {React.Children.map(children, (ch) => (
          <Box display="flex" justifyContent={justifyContent}>
            {ch}
          </Box>
        ))}
      </Stack>
    </Box>
  );
};
