import { useRecoilState } from 'recoil';

import { ParamDropdown } from '@/lib/controls/ParamDropdown';

import { InputRow } from '@/sidebar/ui/InputRow';

import {
  RISK_SOURCE_OPTIONS,
  RISK_TARGET_OPTIONS,
  RiskSourceType,
  RiskTargetType,
  riskSourceTypeState,
  riskTargetTypeState,
} from './state';

export const RiskTypeControl = () => {
  const [riskSource, setRiskSource] = useRecoilState(riskSourceTypeState);
  const [riskTarget, setRiskTarget] = useRecoilState(riskTargetTypeState);

  return (
    <InputRow>
      <ParamDropdown<RiskTargetType>
        title="View"
        value={riskTarget}
        onChange={setRiskTarget}
        options={RISK_TARGET_OPTIONS}
      />
      <ParamDropdown<RiskSourceType>
        title="Affected by"
        value={riskSource}
        onChange={setRiskSource}
        options={RISK_SOURCE_OPTIONS}
      />
    </InputRow>
  );
};
