export const HAZARDS_METADATA = {
  cyclone: {
    label: 'Cyclone',
    dataUnit: 'm/s',
  },
  fluvial: {
    label: 'River Flooding',
    dataUnit: 'm',
  },
  surface: {
    label: 'Surface Flooding',
    dataUnit: 'm',
  },
  coastal: {
    label: 'Coastal Flooding',
    dataUnit: 'm',
  },
};

export const HAZARD_LAYER_NAMES = Object.keys(HAZARDS_METADATA);
