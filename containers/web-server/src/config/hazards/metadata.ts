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
};

export const HAZARDS_MAP_ORDER = ['cyclone', 'fluvial', 'surface', 'coastal'];
export const HAZARDS_UI_ORDER = ['fluvial', 'surface', 'coastal', 'cyclone'];
