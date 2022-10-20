export interface ViewSectionConfig {
  expanded: boolean;
  visible: boolean;
  styles?: string[];
  defaultStyle?: string;
}

export const VIEW_SECTIONS: Record<string, Record<string, ViewSectionConfig>> = {
  hazard: {
    assets: {
      expanded: true,
      visible: true,
      styles: ['type'],
      defaultStyle: 'type',
    },
    hazards: {
      expanded: true,
      visible: true,
    },
    population: {
      expanded: false,
      visible: false,
    },
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
  },
  exposure: {
    assets: {
      expanded: false,
      visible: false,
      styles: ['type'],
      defaultStyle: 'type',
    },
    hazards: {
      expanded: false,
      visible: false,
    },
    population: {
      expanded: true,
      visible: true,
    },
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
  },
  vulnerability: {
    assets: {
      expanded: false,
      visible: false,
      styles: ['type'],
      defaultStyle: 'type',
    },
    hazards: {
      expanded: false,
      visible: false,
    },
    population: {
      expanded: false,
      visible: false,
    },
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
  },
  risk: {
    assets: {
      expanded: false,
      visible: false,

      styles: ['type', 'damages'],
      defaultStyle: 'type',
    },
    hazards: {
      expanded: false,
      visible: true,
    },
    population: {
      expanded: false,
      visible: false,
    },
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
  },
};
