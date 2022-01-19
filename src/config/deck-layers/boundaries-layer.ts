import { MVTLayer } from 'deck.gl';

import { border } from './utils';

export type BoundaryLevel = 'parish' | 'community' | 'subdivision';

function boundaryWidth(level: BoundaryLevel) {
  if (level === 'parish') {
    return 2.5;
  } else if (level === 'community') {
    return 1.5;
  } else if (level === 'subdivision') {
    return 1;
  }
}

export function boundariesLayer(level: BoundaryLevel) {
  return new MVTLayer(
    {
      id: `boundaries-${level}`,
      data: `/vector/data/boundaries_${level}.json`,
      binary: true,
      filled: false,
      refinementStrategy: 'best-available',
      getLineWidth: boundaryWidth(level),
      lineWidthUnits: 'pixels',
    } as any,
    border([190, 190, 190, 255]) as any,
  );
}
