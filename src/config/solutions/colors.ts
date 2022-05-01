import { makeColorConfig } from 'lib/helpers';
import { MarineHabitatType } from './domains';

export const TERRESTRIAL_LANDUSE_COLORS = makeColorConfig({
  Bamboo: '#0f0',
  'Bamboo and Fields': '#0f0',
  'Bamboo and Secondary Forest': '#0f0',
  'Bare Rock': '#333',
  'Bauxite Extraction': '#333',
  'Buildings and other infrastructures': '#900',
  'Closed broadleaved forest (Primary Forest)': '#0f0',
  'Disturbed broadleaved forest (Secondary Forest)': '#0f0',
  'Fields  and Bamboo': '#ff0',
  'Fields and Secondary Forest': '#ff0',
  'Fields or Secondary Forest/Pine Plantation': '#ff0',
  'Fields: Bare Land': '#ff0',
  'Fields: Herbaceous crops, fallow, cultivated vegetables': '#ff0',
  'Fields: Pasture,Human disturbed, grassland': '#ff0',
  'Hardwood Plantation: Euculytus': '##050',
  'Hardwood Plantation: Mahoe': '##050',
  'Hardwood Plantation: Mahogany': '##050',
  'Hardwood Plantation: Mixed': '##050',
  'Herbaceous Wetland': '#0aa',
  'Mangrove Forest': '#044',
  'Open dry forest - Short': '#050',
  'Open dry forest - Tall (Woodland/Savanna)': '#050',
  'Plantation: Tree crops, shrub crops, sugar cane, banana': '#ff0',
  Quarry: '#333',
  'Secondary Forest': '#050',
  'Swamp Forest': '#050',
  'Water Body': '#009',
});

const marineHabitatsConst = [
  {
    value: 'coral',
    label: 'Coral',
  },
  {
    value: 'coral_seagrass',
    label: 'Coral and Seagrass',
  },
  {
    value: 'seagrass',
    label: 'Seagrass',
  },
  {
    value: 'mangrove',
    label: 'Mangrove',
  },
  {
    value: 'coral_mangrove_seagrass',
    label: 'Coral, Mangrove and Seagrass',
  },
  {
    value: 'mangrove_seagrass',
    label: 'Mangrove and Seagrass',
  },
  {
    value: 'coral_mangrove',
    label: 'Coral and Mangrove',
  },
] as const;

export const MARINE_HABITAT_COLORS = makeColorConfig<MarineHabitatType>({
  coral: '#f0f',
  coral_mangrove: '#f00',
  coral_mangrove_seagrass: '#a74',
  coral_seagrass: '#ff0',
  mangrove: '#050',
  mangrove_seagrass: '#080',
  seagrass: '#0f0',
});
