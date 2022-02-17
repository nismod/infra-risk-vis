import { boundariesLayer, BoundaryLevel } from 'config/deck-layers/boundaries-layer';

export function boundariesViewLayer(boundaryLevel: BoundaryLevel) {
  return {
    id: `boundaries_${boundaryLevel}`,
    spatialType: 'vector',
    interactionGroup: 'regions',
    fn: () => boundariesLayer(boundaryLevel),
  };
}
