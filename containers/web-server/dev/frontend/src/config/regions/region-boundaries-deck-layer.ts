import { mvtLayer } from 'lib/deck/layers/base';
import { border } from 'lib/deck/props/style';

import { RegionLevel } from './metadata';

export function regionBoundariesDeckLayer(level: RegionLevel) {
  return mvtLayer(
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
    },
    border([150, 150, 150, 255]),
  );
}
