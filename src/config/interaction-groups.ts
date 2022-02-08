import { makeConfig } from 'lib/helpers';

export const INTERACTION_GROUPS = makeConfig([
  {
    id: 'assets',
    type: 'vector',
    options: {
      pickingRadius: 8,
    },
  },
  {
    id: 'hazards',
    type: 'raster',
  },
  {
    id: 'regions',
    type: 'vector',
  },
]);
