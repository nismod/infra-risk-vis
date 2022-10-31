import { ViewLayer } from '@/lib/data-map/view-layers';

import { RegionLevel } from './metadata';
import { regionBoundariesDeckLayer } from './region-boundaries-deck-layer';

export function regionBoundariesViewLayer(regionLevel: RegionLevel): ViewLayer {
  return {
    id: `boundaries_${regionLevel}`,
    spatialType: 'vector',
    interactionGroup: 'regions',
    params: {
      regionLevel,
    },
    fn: () => regionBoundariesDeckLayer(regionLevel),
  };
}
