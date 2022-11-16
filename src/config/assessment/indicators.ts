import { ValueLabel } from 'lib/controls/params/value-label';

export type IndicatorKey =
  | 'env_ghg'
  | 'env_air_quality'
  | 'env_energy_use'
  | 'env_habitat_disruption'
  | 'env_land'
  | 'econ_passenger'
  | 'econ_freight'
  | 'econ_passenger_occupancy'
  | 'econ_freight_load'
  | 'econ_age'
  | 'econ_road_quality'
  | 'econ_length'
  | 'econ_density'
  | 'econ_border'
  | 'soc_passenger_time'
  | 'soc_passenger_length'
  | 'soc_accidents_death'
  | 'soc_accidents_injury'
  | 'soc_accidents_death_pc'
  | 'soc_accidents_injury_pc'
  | 'soc_noise'
  | 'soc_disease'
  | 'soc_diversity'
  | 'soc_equality'
  | 'soc_inclusivity';

export const INDICATOR_LABELS: ValueLabel<IndicatorKey>[] = [
  { value: 'env_ghg', label: 'GHG emissions' },
  { value: 'env_air_quality', label: 'Air quality' },
  { value: 'env_energy_use', label: 'Energy consumption (non-renewable)' },
  { value: 'env_habitat_disruption', label: 'Habitat and ecosystem disruption' },
  { value: 'env_land', label: 'Land take by transport' },

  { value: 'econ_passenger', label: 'Passenger transport volumes (passenger-km)' },
  { value: 'econ_freight', label: 'Freight transport volumes (tonne-km)' },
  { value: 'econ_passenger_occupancy', label: 'Passenger vehicle occupancy rates' },
  { value: 'econ_freight_load', label: 'Freight transport load factors' },
  { value: 'econ_age', label: 'Average age of vehicle fleet' },
  { value: 'econ_road_quality', label: 'Road quality' },
  { value: 'econ_length', label: 'Length of transport networks' },
  { value: 'econ_density', label: 'Density of transport networks' },
  { value: 'econ_border', label: 'Border restrictions' },

  { value: 'soc_passenger_time', label: 'Average passenger journey time' },
  { value: 'soc_passenger_length', label: 'Average passenger journey length' },
  { value: 'soc_accidents_death', label: 'Total number killed in traffic accidents' },
  { value: 'soc_accidents_injury', label: 'Number of injury traffic accidents' },
  { value: 'soc_accidents_death_pc', label: 'Per capita fatal accident rate' },
  { value: 'soc_accidents_injury_pc', label: 'Per capita injury accident rate' },
  { value: 'soc_noise', label: 'Population affected by traffic noise' },
  { value: 'soc_disease', label: 'Cases of respiratory disease' },
  { value: 'soc_diversity', label: 'Diversity' },
  { value: 'soc_equality', label: 'Equality and fairness' },
  { value: 'soc_inclusivity', label: 'Inclusivity' },
];

export const DEFAULT_WEIGHT = 0.5;
