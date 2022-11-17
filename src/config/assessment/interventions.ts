import { TreeNode } from 'lib/controls/checkbox-tree/tree-node';
import { ValueLabel } from 'lib/controls/params/value-label';
import _ from 'lodash';
import { Effect, ZERO_EFFECT } from './effect';

export type InterventionKey =
  | 'infra_construction'
  | 'infra_maintenance'
  | 'demand_goods'
  | 'demand_travel'
  | 'logistics_planning'
  | 'system_eff'
  | 'fleet_eff'
  | 'fleet_elec'
  | 'road_user_charging'
  | 'custom';

export type InterventionStrength = Record<InterventionKey, number>;
export type InterventionSelection = Record<InterventionKey, boolean>;

export const ZERO_INTERVENTION: InterventionStrength = {
  infra_construction: 0,
  infra_maintenance: 0,
  demand_goods: 0,
  demand_travel: 0,
  logistics_planning: 0,
  system_eff: 0,
  fleet_eff: 0,
  fleet_elec: 0,
  road_user_charging: 0,
  custom: 0,
};

export const NO_INTERVENTIONS: InterventionSelection = {
  infra_construction: false,
  infra_maintenance: false,
  demand_goods: false,
  demand_travel: false,
  logistics_planning: false,
  system_eff: false,
  fleet_eff: false,
  fleet_elec: false,
  road_user_charging: false,
  custom: false,
};

export const INTERVENTION_LABELS: ValueLabel<InterventionKey>[] = [
  { value: 'infra_construction', label: 'Infrastructure construction' },
  { value: 'infra_maintenance', label: 'Infrastructure maintenance' },
  { value: 'demand_goods', label: 'Demand for goods' },
  { value: 'demand_travel', label: 'Demand for travel' },
  { value: 'logistics_planning', label: 'Logistics planning' },
  { value: 'system_eff', label: 'System efficiencies' },
  { value: 'fleet_eff', label: 'Fleet vehicle efficiencies' },
  { value: 'fleet_elec', label: 'Fleet electrification' },
  { value: 'road_user_charging', label: 'Road user charging' },
  { value: 'custom', label: 'Custom intervention' },
];

interface InterventionLabel {
  label: string;
}

export const INTERVENTION_HIERARCHY: TreeNode<InterventionLabel>[] = _.map(
  INTERVENTION_LABELS, ({value, label}: {value: string, label: string}) => ({id: value, label}));

export const INTERVENTION_EFFECTS: Record<InterventionKey, Effect> = {
  infra_construction: {
    ...ZERO_EFFECT,
    env_habitat_disruption: { value: -1 },
    env_land: { value: -1 },
    econ_length: { value: 0.5 },
    econ_density: { value: 0.5 },
    soc_passenger_length: { value: 0.5 },
  },
  infra_maintenance: {
    ...ZERO_EFFECT,
    env_ghg: { value: 0.5 },
    env_energy_use: { value: 0.5 },
    econ_road_quality: { value: 1 },
    soc_passenger_time: { value: 0.5 },
    soc_noise: { value: 0.5 },
  },
  demand_goods: {
    ...ZERO_EFFECT,
    env_ghg: { value: -0.5 },
    env_energy_use: { value: -0.5 },
    econ_freight: { value: 1 },
    econ_freight_load: { value: 1 },
  },
  demand_travel: {
    ...ZERO_EFFECT,
    env_ghg: { value: -0.5 },
    env_energy_use: { value: -0.5 },
    econ_passenger: { value: 1 },
    econ_passenger_occupancy: { value: 1 },
  },
  logistics_planning: {
    ...ZERO_EFFECT,
    econ_freight_load: { value: 0.5 },
    econ_border: { value: 0.5 },
  },
  system_eff: {
    ...ZERO_EFFECT,
    env_ghg: { value: 0.5 },
    env_energy_use: { value: 0.5 },
    econ_freight: { value: 0.5 },
    soc_passenger_time: { value: 0.5 },
  },
  fleet_eff: {
    ...ZERO_EFFECT,
    env_ghg: { value: 0.5 },
    env_energy_use: { value: 0.5 },
  },
  fleet_elec: {
    ...ZERO_EFFECT,
    env_ghg: { value: 1 },
    env_air_quality: { value: 0.5 },
    soc_disease: { value: 0.5 },
  },
  road_user_charging: {
    ...ZERO_EFFECT,
    env_ghg: { value: 0.5 },
    env_energy_use: { value: 0.5 },
    econ_freight: { value: -0.5 },
  },
  custom: ZERO_EFFECT,
};
