import { makeConfig } from 'lib/helpers';
import { InteractionGroupConfig } from 'lib/map/interactions/use-interactions';

export const INTERACTION_GROUPS = makeConfig<InteractionGroupConfig, string>([
  {
    id: 'assets',
    type: 'vector',
    pickingRadius: 8,
    pickMultiple: false,
  },
  {
    id: 'hazards',
    type: 'raster',
    pickMultiple: true,
  },
  {
    id: 'regions',
    type: 'vector',
    pickMultiple: false,
  },
]);
