import { MVTLayer } from 'deck.gl';

import { border } from './utils';

export function boundariesLayer(level: 'parish' | 'community') {
  return new MVTLayer(
    {
      id: `boundaries-${level}`,
      data: `/vector/data/boundaries_${level}.json`,
      binary: true,
      filled: false,
      refinementStrategy: 'best-available',
      getLineWidth: level === 'parish' ? 2 : 1,
      lineWidthUnits: 'pixels',
    } as any,
    border([190, 190, 190, 255]) as any,
  );
}
