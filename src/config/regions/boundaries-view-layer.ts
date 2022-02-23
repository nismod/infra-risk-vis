import { boundariesDeckLayer } from './boundaries-deck-layer';
import { BoundaryLevel } from './metadata';

export function boundariesViewLayer(boundaryLevel: BoundaryLevel) {
  return {
    id: `boundaries_${boundaryLevel}`,
    group: 'regions',
    spatialType: 'vector',
    interactionGroup: 'regions',
    params: {
      boundaryLevel,
    },
    fn: () => boundariesDeckLayer(boundaryLevel),
  };
}
