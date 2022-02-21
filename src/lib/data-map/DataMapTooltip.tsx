import React, { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { hoverPositionState } from './interactions/interaction-state';

export const DataMapTooltip: FC<{}> = ({ children }) => {
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
      {children}
    </div>
  ) : null;
};
