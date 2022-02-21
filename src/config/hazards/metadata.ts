export const HAZARDS_METADATA = {
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
  cyclone: {
    label: 'Cyclone',
    dataUnit: 'km/h', //TODO: check this
  },
};

export const HAZARD_LAYER_NAMES = Object.keys(HAZARDS_METADATA);
