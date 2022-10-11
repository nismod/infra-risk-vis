export const HAZARDS_METADATA = {
  cyclone: {
    label: 'Cyclones',
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
  extreme_heat: {
    label: 'Extreme Heat',
    dataUnit: 'deg',
  },
};

export const HAZARDS_MAP_ORDER = ['cyclone', 'fluvial', 'surface', 'coastal', 'extreme_heat'];
export const HAZARDS_UI_ORDER = ['fluvial', 'surface', 'coastal', 'cyclone', 'extreme_heat'];
