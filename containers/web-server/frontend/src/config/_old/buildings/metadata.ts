import { makeColorConfig, makeConfig } from '@/lib/helpers';

import { AssetMetadata } from '@/config/assets/metadata';

export const BUILDING_COLORS = makeColorConfig({
  commercial: '#f2808c',
  industrial: '#cb97f4',
  institutional: '#808cf2',
  mixed: '#f09e69',
  other: '#bfb4c2',
  recreation: '#95e78b',
  residential: '#f2e680',
  resort: '#f9b2ea',
  unknown: '#dfe4de',
});

export const BUILDING_LAYERS = [
  'buildings_commercial',
  'buildings_industrial',
  'buildings_institutional',
  'buildings_mixed',
  'buildings_other',
  'buildings_recreation',
  'buildings_residential',
  'buildings_resort',
] as const;

export type BuildingLayerType = typeof BUILDING_LAYERS[number];

export const BUILDINGS_METADATA = makeConfig<AssetMetadata, BuildingLayerType>([
  {
    id: 'buildings_commercial',
    type: 'polygon',
    label: 'Buildings (Commercial)',
    color: BUILDING_COLORS.commercial.css,
  },
  {
    id: 'buildings_industrial',
    type: 'polygon',
    label: 'Buildings (Industrial)',
    color: BUILDING_COLORS.industrial.css,
  },
  {
    id: 'buildings_institutional',
    type: 'polygon',
    label: 'Buildings (Institutional)',
    color: BUILDING_COLORS.institutional.css,
  },
  {
    id: 'buildings_mixed',
    type: 'polygon',
    label: 'Buildings (Mixed Use)',
    color: BUILDING_COLORS.mixed.css,
  },
  {
    id: 'buildings_other',
    type: 'polygon',
    label: 'Buildings (Other)',
    color: BUILDING_COLORS.other.css,
  },
  {
    id: 'buildings_recreation',
    type: 'polygon',
    label: 'Buildings (Recreation)',
    color: BUILDING_COLORS.recreation.css,
  },
  {
    id: 'buildings_residential',
    type: 'polygon',
    label: 'Buildings (Residential)',
    color: BUILDING_COLORS.residential.css,
  },
  {
    id: 'buildings_resort',
    type: 'polygon',
    label: 'Buildings (Resort)',
    color: BUILDING_COLORS.resort.css,
  },
]);
