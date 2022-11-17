import _ from 'lodash';
import { atom, selector } from 'recoil';

import { Assessment } from 'config/assessment/assessment';
import { ZERO_EFFECT } from 'config/assessment/effect';
import { DEFAULT_WEIGHT } from 'config/assessment/indicators';
import { INTERVENTION_EFFECTS, INTERVENTION_HIERARCHY, ZERO_INTERVENTION, NO_INTERVENTIONS } from 'config/assessment/interventions';
import { SCENARIO_EFFECTS, ZERO_SCENARIO } from 'config/assessment/scenarios';
import { buildTreeConfig } from 'lib/controls/checkbox-tree/CheckboxTree';

export const currentAssessment = atom<Assessment>({
  key: 'currentAssessment',
  default: {
    description: '',
    notes: '',
    createdAt: new Date(),
    
    interventionSelection: NO_INTERVENTIONS,
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

export const interventionSelection = selector({
  key: 'interventionSelection',
  get: ({ get }) => {
    const assessment = get(currentAssessment);
    return assessment.interventionSelection;
  },
  set: ({ get, set }, newValue) => {
    const assessment = get(currentAssessment);
    set(currentAssessment, {
      ...assessment,
      interventionSelection: newValue,
    });
  },
});

export const interventionTreeConfig = buildTreeConfig(INTERVENTION_HIERARCHY);