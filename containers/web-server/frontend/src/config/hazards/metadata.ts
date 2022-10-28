import { makeOrderingCheck } from '@/lib/helpers';

import { RasterColorMap } from '@/map/legend/RasterLegend';

export const HAZARD_TYPES = ['fluvial', 'coastal', 'cyclone', 'extreme_heat', 'earthquake', 'drought'] as const;

export type HazardType = typeof HAZARD_TYPES[number];

export const HAZARD_COLOR_MAPS: Record<HazardType, RasterColorMap> = {
  fluvial: {
    scheme: 'blues',
    range: [0, 10],
  },
  coastal: {
    scheme: 'greens',
    range: [0, 10],
  },
  cyclone: {
    scheme: 'reds',
    range: [0, 75],
  },
  extreme_heat: {
    scheme: 'reds',
    range: [0, 1],
  },
  earthquake: {
    scheme: 'reds',
    range: [0, 1.4],
  },
  drought: {
    scheme: 'oranges',
    range: [0, 1],
  },
};

export interface HazardMetadata {
  label: string;
  dataUnit: string;
  fractionDigits?: number;
  labelAbbreviations?: Record<string, string>;
  legendAnnotation?: string;
}

export const HAZARDS_METADATA: Record<HazardType, HazardMetadata> = {
  cyclone: {
    label: 'Cyclones',
    dataUnit: 'm/s',
    fractionDigits: 1,
  },
  fluvial: {
    label: 'River Flooding',
    dataUnit: 'm',
    fractionDigits: 1,
  },
  coastal: {
    label: 'Coastal Flooding',
    dataUnit: 'm',
    fractionDigits: 1,
  },
  extreme_heat: {
    label: 'Extreme Heat',
    dataUnit: '',
    legendAnnotation: 'Annual probability of occurrence',
    fractionDigits: 2,
  },
  earthquake: {
    label: 'Seismic Hazard (PGA)',
    dataUnit: 'g',
    fractionDigits: 3,
    labelAbbreviations: {
      PGA: 'Peak Ground Acceleration',
    },
  },
  drought: {
    label: 'Droughts',
    dataUnit: '',
    legendAnnotation: 'Annual probability of occurrence',
    fractionDigits: 2,
  },
};

const hazardOrdering = makeOrderingCheck<HazardType>();

export const HAZARDS_MAP_ORDER = hazardOrdering([
  'earthquake',
  'cyclone',
  'drought',
  'extreme_heat',
  'fluvial',
  'coastal',
]);

export const HAZARDS_UI_ORDER = hazardOrdering([
  'fluvial',
  'coastal',
  'cyclone',
  'drought',
  'extreme_heat',
  'earthquake',
]);
