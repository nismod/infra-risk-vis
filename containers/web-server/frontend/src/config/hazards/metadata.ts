export const HAZARDS_METADATA = {
  cyclone: {
    label: 'Cyclones',
    dataUnit: 'm/s',
  },
  fluvial: {
    label: 'River Flooding',
    dataUnit: 'm',
  },
  coastal: {
    label: 'Coastal Flooding',
    dataUnit: 'm',
  },
  extreme_heat_occurrence: {
    label: 'Extreme Heat (Occurrence)',
    dataUnit: 'leh',
  },
  extreme_heat_exposure: {
    label: 'Extreme Heat (Exposure)',
    dataUnit: 'popn',
  },
  earthquake: {
    label: 'Earthquakes (PGA)',
    dataUnit: 'g',
  },
};

export const HAZARDS_MAP_ORDER = [
  'cyclone',
  'fluvial',
  'coastal',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
  'earthquake',
];
export const HAZARDS_UI_ORDER = [
  'fluvial',
  'coastal',
  'cyclone',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
  'earthquake',
];
