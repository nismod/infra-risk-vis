export interface BoundaryConfig {
  fieldName: string;
  minZoom: number;
  showLabels: boolean;
  labelSingular: string;
  labelPlural: string;
}

export type RegionLevel = 'parish' | 'community' | 'enumeration';

export const REGIONS_METADATA: Record<RegionLevel, BoundaryConfig> = {
  parish: {
    fieldName: 'PARISH',
    minZoom: 9,
    showLabels: true,
    labelSingular: 'Parish',
    labelPlural: 'Parishes',
  },
  community: {
    fieldName: 'COMMUNITY',
    minZoom: 13,
    showLabels: false,
    labelSingular: 'Community',
    labelPlural: 'Communities',
  },
  enumeration: {
    fieldName: 'ED',
    minZoom: 13,
    showLabels: false,
    labelSingular: 'Enumeration District',
    labelPlural: 'Enumeration Districts',
  },
};
