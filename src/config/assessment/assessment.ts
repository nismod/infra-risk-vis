import { Effect, ZERO_EFFECT } from "./effect";
import { InterventionKey, InterventionStrength, InterventionSelection } from "./interventions";
import { ScenarioKey, ScenarioStrength } from "./scenarios";

export interface Assessment {
  description: string;
  notes: string;
  createdAt: Date;

  interventionSelection: InterventionSelection;
  interventionStrength: InterventionStrength;
  defaultInterventionEffects: Record<InterventionKey, Effect>;
  revisedInterventionEffects: Record<InterventionKey, Effect>;

  scenarioStrength: ScenarioStrength;
  defaultScenarioEffects: Record<ScenarioKey, Effect>;
  revisedScenarioEffects: Record<ScenarioKey, Effect>;

  indicatorWeights: Effect;
}

export function weightedSum(unweighted: Effect, weights: Effect, prefix: string = '') {
  let assessed_sum = 0;
  let weighted_sum = 0;
  let count = 0;
  let weight = 0;

  for (let key in unweighted) {
    if (!prefix || key.includes(prefix)) {
      assessed_sum += unweighted[key].value;
      count += 1;
      weighted_sum += unweighted[key].value * weights[key].value;
      weight += weights[key].value;
    }
  }

  const assessed_value = count ? assessed_sum / count : assessed_sum;
  const weighted_value = weight ? weighted_sum / weight : weighted_sum;
  const total_weight = count ? weight / count : count;

  return [assessed_value, weighted_value, total_weight];
}

export function unweightedIndicatorSum(assessment: Assessment): Effect {
  let indicatorsUnweighted: Effect = { ...ZERO_EFFECT };
  for (let key in assessment.revisedScenarioEffects) {
    const effect: Effect = assessment.revisedScenarioEffects[key];
    const strength: number = assessment.scenarioStrength[key];
    for (let indicator in effect) {
      indicatorsUnweighted[indicator] = { value: indicatorsUnweighted[indicator].value + (effect[indicator].value * strength) };
    }
  }
  for (let key in assessment.revisedInterventionEffects) {
    // if the intervention is selected
    if (assessment.interventionSelection[key]){
      // sum up the values of the indicators (in the direction of its strength, +/-1)
      const effect: Effect = assessment.revisedInterventionEffects[key];
      const strength: number = assessment.interventionStrength[key];
      for (let indicator in effect) {
        indicatorsUnweighted[indicator] = { value: indicatorsUnweighted[indicator].value + (effect[indicator].value * strength) };
      }
    }
  }
  return indicatorsUnweighted
}