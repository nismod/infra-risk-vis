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
  {
    value: 'infra_construction',
    label: 'Infrastructure construction',
    description:
      'The construction of new road and rail routes, or expansion of existing routes. This is likely to have an impact on route choice at a local level.',
  },
  {
    value: 'infra_maintenance',
    label: 'Infrastructure maintenance',
    description:
      'Improvements made to existing roads and rail. These are unlikely to have any major effect on transport volumes, but there will be slight capacity improvements.',
  },
  {
    value: 'demand_goods',
    label: 'Demand for goods',
    description:
      'Behaviour change that results in a reduced or increased demand for goods. This could change the requirement for long-distance freight services, with the potential to impact the whole LDT network usage.',
  },
  {
    value: 'demand_travel',
    label: 'Demand for travel',
    description:
      'Any societal behaviour change which affects demand for long-distance travel. This could have a wide-ranging impact on LDT passenger services.',
  },
  {
    value: 'logistics_planning',
    label: 'Logistics planning',
    description:
      'Co-operation between logistics providers and national authorities. The impact that better logistics planning can have on the overall use of the long-distance transport networks is relatively small, but better planning should result in more efficient loading and freight management, together with an easing of restriction at border crossings.',
  },
  {
    value: 'system_eff',
    label: 'System efficiencies',
    description:
      'Technological improvements to vehicle routing and increased use of web-based mobility tools. This will benefit both freight and passenger vehicles on the long-distance transport network, and may help to reduce transport volumes and reduce congestion.',
  },
  {
    value: 'fleet_eff',
    label: 'Fleet vehicle efficiencies',
    description:
      'Improvements to engine efficiencies and use of lighter materials during manufacture. This is likely to reduce emissions and energy consumption very slightly.',
  },
  {
    value: 'fleet_elec',
    label: 'Fleet electrification',
    description:
      'An indication of the number of electric vehicles, particularly freight vehicles using LDT networks. This is likely to be very small in the future, but plans to introduce legislation to phase out fossil-fuelled vehicles are starting to be developed in some sub-Saharan countries.',
  },
  {
    value: 'road_user_charging',
    label: 'Road user charging',
    description:
      'The implementation of long-distance road user charging schemes through tolling interurban toll roads. This could have a significant impact on the long-distance freight movements of specific routes, but as freight companies can pass on any extra costs to customers, there are only likely to be relatively minor impacts on emissions, energy consumption and freight transport volumes.',
  },
  {
    value: 'custom',
    label: 'Custom intervention',
    description: 'An intervention that can be configured to represent any other change to the transport system.',
  },
];

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
