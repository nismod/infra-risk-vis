import { ValueLabel } from '@/lib/controls/params/value-label';

import { RasterColorMap } from '@/map/legend/RasterLegend';

export const BUILDING_DENSITY_TYPES = ['all', 'non_residential'] as const;

export type BuildingDensityType = typeof BUILDING_DENSITY_TYPES[number];

export const BUILDING_DENSITY_TYPE_LABELS: ValueLabel<BuildingDensityType>[] = [
  {
    value: 'all',
    label: 'All',
  },
  {
    value: 'non_residential',
    label: 'Non-residential',
  },
];

export const BUILDING_DENSITY_COLORMAPS: Record<BuildingDensityType, RasterColorMap> = {
  all: {
    scheme: 'orrd',
    range: [0, 500_000],
    rangeTruncated: [false, true],
  },
  non_residential: {
    scheme: 'purples',
    range: [0, 500_000],
    rangeTruncated: [false, true],
  },
};
