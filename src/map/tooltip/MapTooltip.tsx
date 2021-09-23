import React, { FC } from 'react';

interface MapTooltipProps {
  tooltipXY: [number, number];
}

export const MapTooltip: FC<MapTooltipProps> = ({ tooltipXY, children }) => {
  return tooltipXY && React.Children.count(children) ? (
    <div
      style={{
        position: 'absolute',
        zIndex: 1000,
        pointerEvents: 'none',
        left: tooltipXY[0] - 150,
        top: tooltipXY[1] - 5,
      }}
    >
      <div className="tooltip-wrap">
        <div className="tooltip-body">{children}</div>
        <span className="tooltip-triangle"></span>
      </div>
    </div>
  ) : null;
};
