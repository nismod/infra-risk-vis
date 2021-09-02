export const views = {
  overview: {
    layers: [
      'elec_edges_high',
      'elec_edges_low',
      'elec_nodes',
      'rail_edges',
      'rail_nodes',
      'road_edges',
      'bridges',
      'pot_edges',
      'abs_nodes',
    ],
  },
};

export type ViewName = keyof typeof views;
