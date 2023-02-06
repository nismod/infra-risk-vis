import { PropsOf } from '@emotion/react';
import { Box, Stack } from '@mui/material';
import React, { FC } from 'react';

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
    <Box position="absolute" {...effectiveStyle}>
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
