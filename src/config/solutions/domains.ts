const landUseValuesConst = [
  'Bamboo',
  'Bamboo and Fields',
  'Bamboo and Secondary Forest',
  'Bare Rock',
  'Bauxite Extraction',
  'Buildings and other infrastructures',
  'Closed broadleaved forest (Primary Forest)',
  'Disturbed broadleaved forest (Secondary Forest)',
  'Fields  and Bamboo',
  'Fields and Secondary Forest',
  'Fields or Secondary Forest/Pine Plantation',
  'Fields: Bare Land',
  'Fields: Herbaceous crops, fallow, cultivated vegetables',
  'Fields: Pasture,Human disturbed, grassland',
  'Hardwood Plantation: Euculytus',
  'Hardwood Plantation: Mahoe',
  'Hardwood Plantation: Mahogany',
  'Hardwood Plantation: Mixed',
  'Herbaceous Wetland',
  'Mangrove Forest',
  'Open dry forest - Short',
  'Open dry forest - Tall (Woodland/Savanna)',
  'Plantation: Tree crops, shrub crops, sugar cane, banana',
  'Quarry',
  'Secondary Forest',
  'Swamp Forest',
  'Water Body',
] as const;

export const LAND_USE_VALUES = [...landUseValuesConst];

export type LandUseOption = typeof LAND_USE_VALUES[number];

export const TERRESTRIAL_LOCATION_FILTERS = [
  {
    value: 'within_forest_100m',
    label: 'Within 100m of forest',
  },
  {
    value: 'is_protected',
    label: 'Is Protected',
  },
  {
    value: 'is_proposed_protected',
    label: 'Is Proposed Protected',
  },
  {
    value: 'within_major_river_50m',
    label: 'Within 50m of major river',
  },
  {
    value: 'within_large_stream_50m',
    label: 'Within 50m of large stream',
  },
  {
    value: 'within_headwater_stream_50m',
    label: 'Within 50m of headwater stream',
  },
] as const;

export type TerrestrialLocationFilterType = typeof TERRESTRIAL_LOCATION_FILTERS[number]['value'];
