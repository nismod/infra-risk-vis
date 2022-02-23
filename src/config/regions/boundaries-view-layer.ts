import { boundariesDeckLayer, BoundaryLevel } from './boundaries-deck-layer';

export function boundariesViewLayer(boundaryLevel: BoundaryLevel) {
  return {
    id: `boundaries_${boundaryLevel}`,
    group: 'regions',
    spatialType: 'vector',
    interactionGroup: 'regions',
    fn: () => boundariesDeckLayer(boundaryLevel),
  };
}
