import { Box, Paper } from '@mui/material';
import React, { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { hoverPositionState } from 'lib/map/interactions/interaction-state';

export const MapTooltip: FC<{}> = ({ children }) => {
  const tooltipXY = useRecoilValue(hoverPositionState);

  return tooltipXY && React.Children.count(children) ? (
    <div
      style={{
        position: 'absolute',
        zIndex: 1000,
        pointerEvents: 'none',
        left: tooltipXY[0] + 10,
        top: tooltipXY[1],
      }}
    >
      <Paper>
        <Box p={1}>{children}</Box>
      </Paper>
    </div>
  ) : null;
};
