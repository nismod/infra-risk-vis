export interface BoundaryConfig {
  fieldName: string;
  minZoom: number;
  showLabels: boolean;
  labelSingular: string;
  labelPlural: string;
}

export type RegionLevel = 'parish' | 'enumeration';

export const REGIONS_METADATA: Record<RegionLevel, BoundaryConfig> = {
  parish: {
    fieldName: 'PARISH',
    minZoom: 9,
    showLabels: true,
    labelSingular: 'Parish',
    labelPlural: 'Parishes',
  },
  enumeration: {
    fieldName: 'ED',
    minZoom: 13,
    showLabels: false,
    labelSingular: 'Enumeration District',
    labelPlural: 'Enumeration Districts',
  },
};
