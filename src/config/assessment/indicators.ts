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
  {
    value: 'env_ghg',
    label: 'GHG emissions',
    description:
      'Transport-related greenhouse gas (GHG) emissions could be measured in tons of CO₂ per capita. This indicator relates to SDG 9.4. Lower levels of emissions are better.',
  },
  {
    value: 'env_air_quality',
    label: 'Air quality',
    description:
      'Could be measured in terms of PM2.5 air pollution, mean annual exposure. This indicator relates to SDG 3.9. Higher air quality and lower levels of exposure are better.',
  },
  {
    value: 'env_energy_use',
    label: 'Energy consumption (non-renewable)',
    description:
      'Could be measured in terms of non-renewable energy used by mode (ktCO2e). Lower consumption of fossil fuels is better.',
  },
  {
    value: 'env_habitat_disruption',
    label: 'Habitat and ecosystem disruption',
    description:
      'Could be measured in terms of proportion of land area of particular habitat type disrupted by transport infrastructure. Less disruption is better.',
  },
  {
    value: 'env_land',
    label: 'Land take by transport',
    description:
      'Could be measured in terms of proportion of land area required for transport infrastructure. Less land take is better.',
  },

  {
    value: 'econ_passenger',
    label: 'Passenger transport volumes (passenger-km)',
    description:
      'Could be measured in terms of number of passengers. This indicator relates to SDG 9.1. Scoring is context-dependent: less congestion may be an improvement, or greater access may be more important.',
  },
  {
    value: 'econ_freight',
    label: 'Freight transport volumes (tonne-km)',
    description:
      'Could be measured in terms of freight tonnage. This indicator relates to SDG 9.1. Scoring is context-dependent: less congestion may be an improvement, or greater access or trade may be more important.',
  },
  {
    value: 'econ_passenger_occupancy',
    label: 'Passenger vehicle occupancy rates',
    description: 'Could be measured in terms of number of people per vehicle. Higher occupancy is better.',
  },
  {
    value: 'econ_freight_load',
    label: 'Freight transport load factors',
    description:
      'Could be measured in terms of average load factor. More efficient is better (although there could be problems of overloading).',
  },
  {
    value: 'econ_age',
    label: 'Average age of vehicle fleet',
    description:
      'Could be measured in terms of average age in years. Newer vehicles tend to have lower operating emissions, but these have to be traded off against the embodied carbon in vehicle production.',
  },
  {
    value: 'econ_road_quality',
    label: 'Road quality',
    description:
      'Could be measured in terms of an infrastructure quality index. Higher quality is better (although "gold plating" could be an issue).',
  },
  {
    value: 'econ_length',
    label: 'Length of transport networks',
    description:
      'Could be measured in terms of route km of road/rail. Scoring is context-dependent: less congestion is better, but excess capacity is a waste of resources.',
  },
  {
    value: 'econ_density',
    label: 'Density of transport networks',
    description:
      'Could be measured in terms of km of infrastructure per km². Scoring is context-dependent: less congestion and greater access is better, but excess capacity is a waste of resources.',
  },
  {
    value: 'econ_border',
    label: 'Border restrictions',
    description:
      'Could be measured in terms of average delay to freight vehicles at border crossings. Less delay is better.',
  },

  {
    value: 'soc_passenger_time',
    label: 'Average passenger journey time',
    description: 'Could be measured in terms of average journey speed. Faster is better.',
  },
  {
    value: 'soc_passenger_length',
    label: 'Average passenger journey length',
    description:
      'Could be measured in terms of km per average journey. Scoring is context-dependent: unnecessary travel has negative costs, but longer journeys can increase social inclusivity.',
  },
  {
    value: 'soc_accidents_death',
    label: 'Total number killed in traffic accidents',
    description:
      'Could be measured in terms of number of deaths by mode. This indicator relates to SDG 3.6. Fewer deaths are better.',
  },
  {
    value: 'soc_accidents_injury',
    label: 'Number of injury traffic accidents',
    description: 'Could be measured in terms of number of accidents by mode. Fewer accidents are better.',
  },
  {
    value: 'soc_accidents_death_pc',
    label: 'Per capita fatal accident rate',
    description:
      'Could be measured in terms of number of deaths in traffic accidents per capita. Lower rates are better.',
  },
  {
    value: 'soc_accidents_injury_pc',
    label: 'Per capita injury accident rate',
    description:
      'Could be measured in terms of number of injuries in traffic accidents per capita. Lower rates are better.',
  },
  {
    value: 'soc_noise',
    label: 'Population affected by traffic noise',
    description: 'Could be measured in terms of percentage of population affected (by mode). Less exposure is better.',
  },
  {
    value: 'soc_disease',
    label: 'Cases of respiratory disease',
    description:
      'Could be measured in terms of percentage of population with respiratory diseases. Fewer cases are better.',
  },
  {
    value: 'soc_diversity',
    label: 'Diversity',
    description:
      'Could be measured in terms of gender/ethnic split of passengers or labour force. A split which is representative of the local or national population is better.',
  },
  {
    value: 'soc_equality',
    label: 'Equality and fairness',
    description: 'Could be measured in terms of magnitude of gender/ethnic pay gaps. Smaller pay gaps are better.',
  },
  {
    value: 'soc_inclusivity',
    label: 'Inclusivity',
    description:
      'Could be measured in terms of the proportion of population served by intervention who are in the ‘most excluded’ quartile of the national/regional population. Higher inclusivity is better.',
  },
];

export const DEFAULT_WEIGHT = 0.5;
