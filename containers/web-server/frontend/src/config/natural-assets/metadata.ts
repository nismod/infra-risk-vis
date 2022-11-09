import { ValueLabel } from '@/lib/controls/params/value-label';

export const NATURE_RASTER_TYPES = ['biodiversity_intactness', 'forest_landscape_integrity', 'organic_carbon'] as const;

export type NatureRasterType = typeof NATURE_RASTER_TYPES[number];

export const NATURE_RASTER_VALUE_LABELS: ValueLabel<NatureRasterType>[] = [
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
    label: 'Soil Organic Carbon',
  },
];
