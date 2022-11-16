import { ValueLabel } from 'lib/controls/params/value-label';
import { Effect, ZERO_EFFECT } from './effect';

export type ScenarioKey = 'population' | 'economic' | 'energy-cost';

export type ScenarioStrength = Record<ScenarioKey, number>;

export const ZERO_SCENARIO: ScenarioStrength = {
  population: 0,
  economic: 0,
  'energy-cost': 0,
};

export const SCENARIO_LABELS: ValueLabel<ScenarioKey>[] = [
  { value: 'population', label: 'Population' },
  { value: 'economic', label: 'Economic' },
  { value: 'energy-cost', label: 'Energy Costs' },
];

export const SCENARIO_EFFECTS: Record<ScenarioKey, Effect> = {
  population: {
    ...ZERO_EFFECT,
    env_ghg: { value: -1 },
    env_energy_use: { value: -1 },
    econ_passenger: { value: 1 },
    econ_freight: { value: 1 },
    soc_accidents_death: { value: -0.5 },
    soc_accidents_injury: { value: -0.5 },
  },
  economic: {
    ...ZERO_EFFECT,
    env_ghg: { value: -1 },
    env_energy_use: { value: -1 },
    econ_passenger: { value: 1 },
    econ_freight: { value: 1 },
    econ_age: { value: 0.5 },
  },
  'energy-cost': {
    ...ZERO_EFFECT,
    env_ghg: { value: 1 },
    env_energy_use: { value: 1 },
    econ_passenger: { value: -1 },
    econ_freight: { value: -1 },
  },
};
