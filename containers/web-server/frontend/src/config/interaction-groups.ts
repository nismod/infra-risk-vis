import { InteractionGroupConfig } from '@/lib/data-map/interactions/use-interactions';
import { makeConfig } from '@/lib/helpers';

export const INTERACTION_GROUPS = makeConfig<InteractionGroupConfig, string>([
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
    id: 'population',
    type: 'raster',
    pickMultiple: false,
  },
  {
    id: 'hdi',
    type: 'vector',
    usesAutoHighlight: true,
    pickMultiple: false,
  },
  {
    id: 'regions',
    type: 'vector',
    pickingRadius: 8,
    pickMultiple: false,
  },
  {
    id: 'solutions',
    type: 'vector',
    pickingRadius: 8,
    usesAutoHighlight: true,
    pickMultiple: false,
  },
  {
    id: 'drought',
    type: 'vector',
    pickingRadius: 8,
    usesAutoHighlight: true,
    pickMultiple: false,
  },
  {
    id: 'wdpa',
    type: 'vector',
    pickingRadius: 0,
    usesAutoHighlight: true,
    pickMultiple: true,
  },
  {
    id: 'buildings_ghsl',
    type: 'raster',
    pickMultiple: false,
  },
]);
