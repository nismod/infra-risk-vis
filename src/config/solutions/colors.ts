import { makeColorConfig } from 'lib/helpers';
import { MarineHabitatType } from './domains';

export const TERRESTRIAL_LANDUSE_COLORS = makeColorConfig({
  Bamboo: '#8cae26',
  'Bamboo and Fields': '#b2c73d',
  'Bamboo and Secondary Forest': '#6caa10',
  'Bare Rock': '#dcdcdc',
  'Bauxite Extraction': '#9a5fa9',
  'Buildings and other infrastructures': '#c43c39',
  'Closed broadleaved forest (Primary Forest)': '#117733',
  'Disturbed broadleaved forest (Secondary Forest)': '#44aa66',
  'Fields  and Bamboo': '#b2c73d',
  'Fields and Secondary Forest': '#e6cc4d',
  'Fields or Secondary Forest/Pine Plantation': '#c9a063',
  'Fields: Bare Land': '#ffee88',
  'Fields: Herbaceous crops, fallow, cultivated vegetables': '#ffffa8',
  'Fields: Pasture,Human disturbed, grassland': '#e6e64d',
  'Hardwood Plantation: Euculytus': '#aa1636',
  'Hardwood Plantation: Mahoe': '#ac495d',
  'Hardwood Plantation: Mahogany': '#862a3b',
  'Hardwood Plantation: Mixed': '#8a3d4e',
  'Herbaceous Wetland': '#6beabc',
  'Mangrove Forest': '#80f2e6',
  'Open dry forest - Short': '#49df51',
  'Open dry forest - Tall (Woodland/Savanna)': '#82e47f',
  'Plantation: Tree crops, shrub crops, sugar cane, banana': '#f2a64d',
  Quarry: '#bf8acb',
  'Secondary Forest': '#73ba6b',
  'Swamp Forest': '#37d47e',
  'Water Body': '#6bdbe9',
});

export const MARINE_HABITAT_COLORS = makeColorConfig<MarineHabitatType>({
  coral: '#f9b2ea',
  coral_mangrove: '#808cf2',
  coral_mangrove_seagrass: '#152513',
  coral_seagrass: '#f2808c',
  mangrove: '#80f2e6',
  mangrove_seagrass: '#95e78b',
  seagrass: '#f2e680',
  other: 'rgba(128,128,128,0.2)',
});
