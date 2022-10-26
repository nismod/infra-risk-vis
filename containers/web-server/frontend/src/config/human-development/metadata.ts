import { ValueLabel } from '@/lib/controls/params/value-label';

export const HDI_REGION_LEVELS = ['countries', 'regions'] as const;

export type HdiRegionLevel = typeof HDI_REGION_LEVELS[number];

export const HDI_REGION_LEVEL_LABELS: ValueLabel<HdiRegionLevel>[] = [
  {
    value: 'countries',
    label: 'Countries',
  },
  {
    value: 'regions',
    label: 'Regions',
  },
];

export const HDI_REGION_LEVEL_METADATA = {
  countries: {
    nameField: 'country',
  },
  regions: {
    nameField: 'region',
  },
};

export const HDI_VARIABLES = ['subnational_hdi', 'health_index', 'educational_index', 'income_index'] as const;

export type HdiVariableType = typeof HDI_VARIABLES[number];

export const HDI_VARIABLE_LABELS: ValueLabel<HdiVariableType>[] = [
  {
    value: 'subnational_hdi',
    label: 'Human Development Index',
  },
  {
    value: 'health_index',
    label: 'Health Index',
  },
  {
    value: 'educational_index',
    label: 'Educational Index',
  },
  {
    value: 'income_index',
    label: 'Income Index',
  },
];
