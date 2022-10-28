export const HAZARD_COLOR_MAPS = {
  fluvial: {
    scheme: 'blues',
    range: [0, 10],
  },
  coastal: {
    scheme: 'greens',
    range: [0, 10],
  },
  surface: {
    scheme: 'purples',
    range: [0, 10],
  },
  cyclone: {
    scheme: 'reds',
    range: [0, 75],
  },
  extreme_heat_exposure: {
    scheme: 'reds',
    range: [0, 250],
  },
  extreme_heat_occurrence: {
    scheme: 'reds',
    range: [0, 1],
  },
  earthquake: {
    scheme: 'reds',
    range: [0, 1.4],
  },
  drought: {
    scheme: 'oranges',
    range: [0, 100],
  },
};

export const HAZARDS_METADATA = {
  cyclone: {
    label: 'Cyclones',
    dataUnit: 'm/s',
    fractionDigits: 1,
  },
  fluvial: {
    label: 'River Flooding',
    dataUnit: 'm',
    fractionDigits: 1,
  },
  coastal: {
    label: 'Coastal Flooding',
    dataUnit: 'm',
    fractionDigits: 1,
  },
  extreme_heat_occurrence: {
    label: 'Extreme Heat (Occurrence)',
    dataUnit: 'leh',
    fractionDigits: 1,
  },
  extreme_heat_exposure: {
    label: 'Extreme Heat (Exposure)',
    dataUnit: 'popn',
    fractionDigits: 0,
  },
  earthquake: {
    label: 'Seismic Hazard (PGA)',
    dataUnit: 'g',
    fractionDigits: 3,
    labelAbbreviations: {
      PGA: 'Peak Ground Acceleration',
    },
  },
  drought: {
    label: 'Droughts',
    dataUnit: '',
    // fractionDigits:0
  },
};

export const HAZARDS_MAP_ORDER = [
  'earthquake',
  'cyclone',
  'drought',
  'fluvial',
  'coastal',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
];
export const HAZARDS_UI_ORDER = [
  'fluvial',
  'coastal',
  'cyclone',
  'drought',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
  'earthquake',
];
