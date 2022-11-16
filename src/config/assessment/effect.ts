import { IndicatorKey } from './indicators';

export interface AnnotatedValue {
  value: number;
  notes?: string;
}

export type Effect = Record<IndicatorKey, AnnotatedValue>;

export const ZERO_EFFECT = {
  env_ghg: { value: 0 },
  env_air_quality: { value: 0 },
  env_energy_use: { value: 0 },
  env_habitat_disruption: { value: 0 },
  env_land: { value: 0 },
  econ_passenger: { value: 0 },
  econ_freight: { value: 0 },
  econ_passenger_occupancy: { value: 0 },
  econ_freight_load: { value: 0 },
  econ_age: { value: 0 },
  econ_road_quality: { value: 0 },
  econ_length: { value: 0 },
  econ_density: { value: 0 },
  econ_border: { value: 0 },
  soc_passenger_time: { value: 0 },
  soc_passenger_length: { value: 0 },
  soc_accidents_death: { value: 0 },
  soc_accidents_injury: { value: 0 },
  soc_accidents_death_pc: { value: 0 },
  soc_accidents_injury_pc: { value: 0 },
  soc_noise: { value: 0 },
  soc_disease: { value: 0 },
  soc_diversity: { value: 0 },
  soc_equality: { value: 0 },
  soc_inclusivity: { value: 0 },
};
