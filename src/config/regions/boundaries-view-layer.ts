import { regionBoundariesDeckLayer } from './region-boundaries-deck-layer';
import { RegionLevel } from './metadata';

export function regionBoundariesViewLayer(regionLevel: RegionLevel) {
  return {
    id: `region_boundaries_${regionLevel}`,
    group: 'regions',
    spatialType: 'vector',
    interactionGroup: 'regions',
    params: {
      regionLevel,
    },
    fn: () => regionBoundariesDeckLayer(regionLevel),
  };
}
