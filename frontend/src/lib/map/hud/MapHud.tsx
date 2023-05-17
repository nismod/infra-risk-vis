import { Box } from '@mui/system';
import { FC } from 'react';

export const MapHud: FC<{ top?: number; right?: number; bottom?: number; left?: number }> = ({
  children,
  top = 0,
  right = 0,
  bottom = 0,
  left = 0,
}) => (
  <Box
    position="absolute"
    {...{ top, right, bottom, left }}
    zIndex={1000}
    sx={{ pointerEvents: 'none' }}
  >
    {children}
  </Box>
);
