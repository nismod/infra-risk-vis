import { ValueLabel } from '@/lib/controls/params/value-label';

export const REGIONAL_EXPOSURE_VARIABLES = [
  'pop_exposed_seismic_threshold0.1g',
  'pop_exposed_seismic_threshold0.2g',
  'pop_exposed_river_historical_WATCH_1980_thresholdNone',
  'pop_exposed_river_rcp4p5_MIROC-ESM-CHEM_2050_thresholdNone',
] as const;

export type RegionalExposureVariableType = typeof REGIONAL_EXPOSURE_VARIABLES[number];

export const REGIONAL_EXPOSURE_VARIABLE_LABELS: ValueLabel<RegionalExposureVariableType>[] = [
  {
    value: 'pop_exposed_seismic_threshold0.1g',
    label: 'Seismic hazard >0.1g',
  },
  {
    value: 'pop_exposed_seismic_threshold0.2g',
    label: 'Seismic hazard >0.2g',
  },
  {
    value: 'pop_exposed_river_historical_WATCH_1980_thresholdNone',
    label: 'River flooding (baseline, 100yr)'
  },
  {
    value: 'pop_exposed_river_rcp4p5_MIROC-ESM-CHEM_2050_thresholdNone',
    label: 'River flooding (RCP4.5, 2050, 100yr)'
  }
];
