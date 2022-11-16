import _ from 'lodash';
import { atom, selector } from 'recoil';

import { Assessment } from 'config/assessment/assessment';
import { ZERO_EFFECT } from 'config/assessment/effect';
import { DEFAULT_WEIGHT } from 'config/assessment/indicators';
import { INTERVENTION_EFFECTS, ZERO_INTERVENTION } from 'config/assessment/interventions';
import { SCENARIO_EFFECTS, ZERO_SCENARIO } from 'config/assessment/scenarios';

export const currentAssessment = atom<Assessment>({
  key: 'currentAssessment',
  default: {
    description: '',
    notes: '',
    createdAt: new Date(),

    interventionStrength: ZERO_INTERVENTION,
    defaultInterventionEffects: INTERVENTION_EFFECTS,
    revisedInterventionEffects: INTERVENTION_EFFECTS,

    scenarioStrength: ZERO_SCENARIO,
    defaultScenarioEffects: SCENARIO_EFFECTS,
    revisedScenarioEffects: SCENARIO_EFFECTS,

    indicatorWeights: _.mapValues(ZERO_EFFECT, () => ({ value: DEFAULT_WEIGHT })),
  },
});

export const indicatorWeights = selector({
  key: 'indicatorWeights',
  get: ({ get }) => {
    const assessment = get(currentAssessment);
    return assessment.indicatorWeights;
  },
  set: ({ get, set }, newValue) => {
    const assessment = get(currentAssessment);
    set(currentAssessment, {
      ...assessment,
      indicatorWeights: newValue,
    });
  },
});
