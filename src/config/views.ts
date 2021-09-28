import { makeConfig } from '../helpers';

export const VIEWS = makeConfig([
  {
    id: 'overview',
    layers: [
      'hazard_fluvial_20',
      'hazard_fluvial_50',
      'hazard_fluvial_100',
      'hazard_fluvial_200',
      'hazard_fluvial_500',
      'hazard_fluvial_1500',
      'hazard_coastal_1',
      'hazard_coastal_2',
      'hazard_coastal_5',
      'hazard_coastal_10',
      'hazard_coastal_50',
      'hazard_coastal_100',
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
  {
    id: 'hazards',
    layers: [
      'hazard_fluvial_20',
      'hazard_fluvial_50',
      'hazard_fluvial_100',
      'hazard_fluvial_200',
      'hazard_fluvial_500',
      'hazard_fluvial_1500',
      'hazard_coastal_1',
      'hazard_coastal_2',
      'hazard_coastal_5',
      'hazard_coastal_10',
      'hazard_coastal_50',
      'hazard_coastal_100',
    ],
  },
]);

export type ViewName = keyof typeof VIEWS;
