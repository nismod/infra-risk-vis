import { InteractionGroupConfig } from 'lib/data-map/interactions/use-interactions';

export const INTERACTION_GROUPS: InteractionGroupConfig[] = [
  {
    id: 'assets',
    type: 'vector',
    pickingRadius: 8,
    pickMultiple: false,
    usesAutoHighlight: true,
  },
  {
    id: 'hazards',
    type: 'raster',
    pickMultiple: true,
  },
  {
    id: 'regions',
    type: 'vector',
    pickingRadius: 8,
    pickMultiple: false,
  },
];
