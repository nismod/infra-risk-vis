import { TreeNode } from 'lib/controls/checkbox-tree/tree-node';

interface LanduseTypeData {
  id: string;
  label: string;
}

export const LANDUSE_HIERARCHY: TreeNode<LanduseTypeData>[] = [
  {
    id: 'some-group',
    label: 'Some Group',
    children: [
      { id: 'Bamboo', label: 'Bamboo' },
      { id: 'Bamboo and Fields', label: 'Bamboo and Fields' },
      { id: 'Bamboo and Secondary Forest', label: 'Bamboo and Secondary Forest' },
      { id: 'Bare Rock', label: 'Bare Rock' },
      { id: 'Bauxite Extraction', label: 'Bauxite Extraction' },
      { id: 'Buildings and other infrastructures', label: 'Buildings and other infrastructures' },
      { id: 'Closed broadleaved forest (Primary Forest)', label: 'Closed broadleaved forest (Primary Forest)' },
      {
        id: 'Disturbed broadleaved forest (Secondary Forest)',
        label: 'Disturbed broadleaved forest (Secondary Forest)',
      },
      { id: 'Fields  and Bamboo', label: 'Fields  and Bamboo' },
      { id: 'Fields and Secondary Forest', label: 'Fields and Secondary Forest' },
      { id: 'Fields or Secondary Forest/Pine Plantation', label: 'Fields or Secondary Forest/Pine Plantation' },
      { id: 'Fields: Bare Land', label: 'Fields: Bare Land' },
      {
        id: 'Fields: Herbaceous crops, fallow, cultivated vegetables',
        label: 'Fields: Herbaceous crops, fallow, cultivated vegetables',
      },
    ],
  },
  {
    id: 'some-other-group',
    label: 'Some Other Group',
    children: [
      { id: 'Fields: Pasture,Human disturbed, grassland', label: 'Fields: Pasture,Human disturbed, grassland' },
      { id: 'Hardwood Plantation: Euculytus', label: 'Hardwood Plantation: Euculytus' },
      { id: 'Hardwood Plantation: Mahoe', label: 'Hardwood Plantation: Mahoe' },
      { id: 'Hardwood Plantation: Mahogany', label: 'Hardwood Plantation: Mahogany' },
      { id: 'Hardwood Plantation: Mixed', label: 'Hardwood Plantation: Mixed' },
      { id: 'Herbaceous Wetland', label: 'Herbaceous Wetland' },
      { id: 'Mangrove Forest', label: 'Mangrove Forest' },
      { id: 'Open dry forest - Short', label: 'Open dry forest - Short' },
      { id: 'Open dry forest - Tall (Woodland/Savanna)', label: 'Open dry forest - Tall (Woodland/Savanna)' },
      {
        id: 'Plantation: Tree crops, shrub crops, sugar cane, banana',
        label: 'Plantation: Tree crops, shrub crops, sugar cane, banana',
      },
      { id: 'Quarry', label: 'Quarry' },
      { id: 'Secondary Forest', label: 'Secondary Forest' },
      { id: 'Swamp Forest', label: 'Swamp Forest' },
      { id: 'Water Body', label: 'Water Body' },
    ],
  },
];
