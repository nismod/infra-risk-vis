import { regionBoundariesDeckLayer } from './region-boundaries-deck-layer';
import { RegionLevel } from './metadata';
import { ViewLayer } from 'lib/data-map/view-layers';

export function regionBoundariesViewLayer(regionLevel: RegionLevel): ViewLayer {
  return {
    id: `boundaries_${regionLevel}`,
    group: 'regions',
    spatialType: 'vector',
    interactionGroup: 'regions',
    params: {
      regionLevel,
    },
    fn: () => regionBoundariesDeckLayer(regionLevel),
  };
}
