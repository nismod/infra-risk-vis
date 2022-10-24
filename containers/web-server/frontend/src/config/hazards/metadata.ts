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
};

export const HAZARDS_MAP_ORDER = [
  'earthquake',
  'cyclone',
  'fluvial',
  'coastal',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
];
export const HAZARDS_UI_ORDER = [
  'fluvial',
  'coastal',
  'cyclone',
  'extreme_heat_occurrence',
  'extreme_heat_exposure',
  'earthquake',
];
