import { TreeNode } from '@/lib/controls/checkbox-tree/tree-node';

interface LanduseTypeData {
  id: string;
  label: string;
}

export const LANDUSE_HIERARCHY: TreeNode<LanduseTypeData>[] = [
  {
    id: 'agriculture',
    label: 'Agriculture',
    children: [
      { id: 'Fields: Bare Land', label: 'Fields: Bare Land' },
      {
        id: 'Fields: Herbaceous crops, fallow, cultivated vegetables',
        label: 'Fields: Herbaceous crops, fallow, cultivated vegetables',
      },
      { id: 'Fields: Pasture,Human disturbed, grassland', label: 'Fields: Pasture, human-disturbed, grassland' },
    ],
  },
  { id: 'Bare Rock', label: 'Bare Rock' },
  {
    id: 'bamboo',
    label: 'Bamboo and Mixed',
    children: [
      { id: 'Bamboo', label: 'Bamboo' },
      { id: 'Bamboo and Fields', label: 'Bamboo and Fields' },
      { id: 'Bamboo and Secondary Forest', label: 'Bamboo and Secondary Forest' },
      { id: 'Fields  and Bamboo', label: 'Fields and Bamboo' },
    ],
  },
  { id: 'Buildings and other infrastructures', label: 'Built-up areas' },
  {
    id: 'forest',
    label: 'Forest',
    children: [
      { id: 'Closed broadleaved forest (Primary Forest)', label: 'Closed broadleaved forest (Primary Forest)' },
      {
        id: 'Disturbed broadleaved forest (Secondary Forest)',
        label: 'Disturbed broadleaved forest (Secondary Forest)',
      },
      { id: 'Secondary Forest', label: 'Secondary Forest' },
      { id: 'Swamp Forest', label: 'Swamp Forest' },
      { id: 'Mangrove Forest', label: 'Mangrove Forest' },
    ],
  },
  {
    id: 'mining',
    label: 'Mining',
    children: [
      { id: 'Bauxite Extraction', label: 'Bauxite Extraction' },
      { id: 'Quarry', label: 'Quarry' },
    ],
  },
  {
    id: 'mixed',
    label: 'Mixed Use',
    children: [
      { id: 'Fields and Secondary Forest', label: 'Fields and Secondary Forest' },
      { id: 'Fields or Secondary Forest/Pine Plantation', label: 'Fields or Secondary Forest/Pine Plantation' },
    ],
  },
  {
    id: 'open-dry-forest',
    label: 'Open dry forest',
    children: [
      { id: 'Open dry forest - Short', label: 'Open dry forest - Short' },
      { id: 'Open dry forest - Tall (Woodland/Savanna)', label: 'Open dry forest - Tall (Woodland/Savanna)' },
    ],
  },
  {
    id: 'plantation',
    label: 'Plantation',
    children: [
      { id: 'Hardwood Plantation: Euculytus', label: 'Hardwood Plantation: Euculytus' },
      { id: 'Hardwood Plantation: Mahoe', label: 'Hardwood Plantation: Mahoe' },
      { id: 'Hardwood Plantation: Mahogany', label: 'Hardwood Plantation: Mahogany' },
      { id: 'Hardwood Plantation: Mixed', label: 'Hardwood Plantation: Mixed' },
      {
        id: 'Plantation: Tree crops, shrub crops, sugar cane, banana',
        label: 'Plantation: Tree crops, shrub crops, sugar cane, banana',
      },
    ],
  },
  { id: 'Water Body', label: 'Water Body' },
  {
    id: 'wetland',
    label: 'Wetland',
    children: [{ id: 'Herbaceous Wetland', label: 'Herbaceous Wetland' }],
  },
];
