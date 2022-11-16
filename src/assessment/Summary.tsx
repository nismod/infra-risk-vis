import { weightedSum } from "config/assessment/assessment";
import { Effect } from "config/assessment/effect";
import { useRecoilState } from "recoil";
import { indicatorWeights } from "state/assessment";
import { QualitativeText } from "./QualitativeText";

export const Summary = ({
  unweighted,
  overall_assessed,
  overall_weighted,
}: {
  unweighted: Effect;
  overall_assessed: Number;
  overall_weighted: Number;
}) => {
  const [currentWeights, _setWeights] = useRecoilState(indicatorWeights);

  // Force our way around type-checking - recoil returns the expected Effect
  // @ts-ignore
  let weights: Effect = { ...currentWeights };

  const [_assessed_env, weighted_env, _total_weight_env] = weightedSum(unweighted, weights, 'env');
  const [_assessed_soc, weighted_soc, _total_weight_soc] = weightedSum(unweighted, weights, 'soc');
  const [_assessed_econ, weighted_econ, _total_weight_econ] = weightedSum(unweighted, weights, 'econ');

  return (
    <>
      <h2>Summary</h2>
      <p>Overall, the proposed interventions are expected to have:</p>
      <ul>
        <li>
          a <QualitativeText value={weighted_env} /> effect on environmental sustainability
        </li>
        <li>
          a <QualitativeText value={weighted_soc} /> effect on social sustainability
        </li>
        <li>
          a <QualitativeText value={weighted_econ} /> effect on economic sustainability
        </li>
      </ul>
      <p>
        Given the weights assigned, this could be consider a <QualitativeText value={overall_weighted} />{' '}
        effect overall.
      </p>
    </>
  );
};