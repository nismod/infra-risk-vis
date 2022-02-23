export interface BoundaryConfig {
  fieldName: string;
  minZoom: number;
  showLabels: boolean;
  label: string;
}

export type BoundaryLevel = 'parish' | 'community' | 'enumeration';

export const REGIONS_METADATA: Record<BoundaryLevel, BoundaryConfig> = {
  parish: {
    fieldName: 'PARISH',
    minZoom: 9,
    showLabels: true,
    label: 'Parish',
  },
  community: {
    fieldName: 'COMMUNITY',
    minZoom: 13,
    showLabels: false,
    label: 'Community',
  },
  enumeration: {
    fieldName: 'ED',
    minZoom: 13,
    showLabels: false,
    label: 'Enumeration District',
  },
};
