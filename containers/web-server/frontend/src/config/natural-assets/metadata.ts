import { ValueLabel } from '@/lib/controls/params/value-label';

export const NATURAL_ASSET_TYPES = ['biodiversity_intactness', 'forest_landscape_integrity', 'organic_carbon'] as const;

export type NaturalAssetType = typeof NATURAL_ASSET_TYPES[number];

export const NATURAL_ASSET_VALUE_LABELS: ValueLabel<NaturalAssetType>[] = [
  {
    value: 'biodiversity_intactness',
    label: 'Biodiversity Intactness',
  },
  {
    value: 'forest_landscape_integrity',
    label: 'Forest Landscape Integrity',
  },
  {
    value: 'organic_carbon',
    label: 'Organic Carbon',
  },
];
