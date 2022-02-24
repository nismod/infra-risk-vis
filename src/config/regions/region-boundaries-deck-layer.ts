import { MVTLayer } from 'deck.gl';

import { border } from 'lib/deck-layers/utils';

import { RegionLevel } from './metadata';

export function regionBoundariesDeckLayer(level: RegionLevel) {
  return new MVTLayer(
    {
      id: `boundaries_${level}`,
      data: `/vector/data/regions_${level}.json`,
      binary: true,
      filled: true,
      getFillColor: [255, 255, 255, 0],
      pickable: true,
      stroked: true,
      refinementStrategy: 'best-available',
      lineWidthUnits: 'pixels',
    } as any,
    border([150, 150, 150, 255]) as any,
  );
}
