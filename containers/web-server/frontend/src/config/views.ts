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
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
    regions: {
      expanded: false,
      visible: false,

      styles: ['boundaries', 'population'],
      defaultStyle: 'boundaries',
    },
  },
  exposure: {
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
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
    regions: {
      expanded: false,
      visible: false,

      styles: ['boundaries', 'population'],
      defaultStyle: 'boundaries',
    },
  },
  vulnerability: {
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
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
    regions: {
      expanded: false,
      visible: false,

      styles: ['boundaries', 'population'],
      defaultStyle: 'boundaries',
    },
  },
  risk: {
    assets: {
      expanded: true,
      visible: true,

      styles: ['type', 'damages'],
      defaultStyle: 'damages',
    },
    hazards: {
      expanded: false,
      visible: true,
    },
    buildings: {
      expanded: false,
      visible: false,

      styles: ['type'],
      defaultStyle: 'type',
    },
    regions: {
      expanded: false,
      visible: false,

      styles: ['boundaries', 'population'],
      defaultStyle: 'boundaries',
    },
    terrestrial: {
      expanded: false,
      visible: false,
      styles: ['landuse', 'slope', 'elevation'],
      defaultStyle: 'landuse',
    },
    marine: {
      expanded: false,
      visible: false,
      styles: ['habitat'],
      defaultStyle: 'habitat',
    },
  },
};
