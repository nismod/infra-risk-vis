import { makeConfig } from '@/lib/helpers';

export const TERRESTRIAL_STYLES = makeConfig([
  { id: 'landuse', label: 'Land Use' },
  { id: 'slope', label: 'Slope' },
  { id: 'elevation', label: 'Elevation' },
]);

export const MARINE_STYLES = makeConfig([{ id: 'habitat', label: 'Habitat' }]);
