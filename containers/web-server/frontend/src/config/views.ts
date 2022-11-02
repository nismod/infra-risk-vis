export interface ViewSectionConfig {
  expanded: boolean;
  visible: boolean;
  styles?: string[];
  defaultStyle?: string;
}

export const VIEW_SECTIONS: Record<string, Record<string, ViewSectionConfig>> = {
  hazard: {
    assets: {
      expanded: false,
      visible: false,
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
    },
    industry: {
      expanded: false,
      visible: false,
    },
    healthcare: {
      expanded: false,
      visible: false,
    },
    'natural-assets': {
      expanded: false,
      visible: false,
    },
    'human-vulnerability': {
      expanded: false,
      visible: false,
    },
    'travel-time': {
      expanded: false,
      visible: false,
    },
    'nature-vulnerability': {
      expanded: false,
      visible: false,
    },
  },
  exposure: {
    assets: {
      expanded: true,
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
    },
    industry: {
      expanded: true,
      visible: false,
    },
    healthcare: {
      expanded: true,
      visible: false,
    },
    'natural-assets': {
      expanded: true,
      visible: false,
    },
    'human-vulnerability': {
      expanded: false,
      visible: false,
    },
    'travel-time': {
      expanded: false,
      visible: false,
    },
    'nature-vulnerability': {
      expanded: false,
      visible: false,
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
    },
    industry: {
      expanded: false,
      visible: false,
    },
    healthcare: {
      expanded: false,
      visible: false,
    },
    'natural-assets': {
      expanded: false,
      visible: false,
    },
    'human-vulnerability': {
      expanded: true,
      visible: true,
    },
    'travel-time': {
      expanded: true,
      visible: false,
    },
    'nature-vulnerability': {
      expanded: true,
      visible: false,
    },
  },
  risk: {
    assets: {
      expanded: false,
      visible: false,

      styles: ['type'],
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
    },
    industry: {
      expanded: false,
      visible: false,
    },
    healthcare: {
      expanded: false,
      visible: false,
    },
    'natural-assets': {
      expanded: false,
      visible: false,
    },
    'human-vulnerability': {
      expanded: false,
      visible: false,
    },
    'travel-time': {
      expanded: false,
      visible: false,
    },
    'nature-vulnerability': {
      expanded: false,
      visible: false,
    },
  },
};
