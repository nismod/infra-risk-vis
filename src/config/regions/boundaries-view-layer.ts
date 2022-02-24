import { regionBoundariesDeckLayer } from './region-boundaries-deck-layer';
import { RegionLevel } from './metadata';

export function regionBoundariesViewLayer(boundaryLevel: RegionLevel) {
  return {
    id: `boundaries_${boundaryLevel}`,
    group: 'regions',
    spatialType: 'vector',
    interactionGroup: 'regions',
    params: {
      boundaryLevel,
    },
    fn: () => regionBoundariesDeckLayer(boundaryLevel),
  };
}
