export interface BoundaryConfig {
  fieldName: string;
  minZoom: number;
  showLabels: boolean;
  labelSingular: string;
  labelPlural: string;
}

export type RegionLevel = 'level0' | 'level1';

export const REGIONS_METADATA: Record<RegionLevel, BoundaryConfig> = {
  level0: {
    fieldName: 'NAME_0',
    minZoom: 3,
    showLabels: true,
    labelSingular: 'Country',
    labelPlural: 'Countries',
  },
  level1: {
    fieldName: 'NAME_1',
    minZoom: 3,
    showLabels: true,
    labelSingular: 'Admin-1 Region',
    labelPlural: 'Admin-1 Regions',
  },
};
