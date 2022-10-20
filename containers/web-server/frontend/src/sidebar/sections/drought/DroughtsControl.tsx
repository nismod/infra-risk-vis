import { Checkbox, Divider, FormControlLabel } from '@mui/material';
import { FC } from 'react';
import { useRecoilState } from 'recoil';

import { ParamDropdown } from '@/lib/controls/ParamDropdown';

import {
  DROUGHT_OPTIONS_VARIABLE_LABELS,
  DROUGHT_RISK_VARIABLE_LABELS,
  DroughtOptionsVariableType,
  DroughtRiskVariableType,
} from '@/config/drought/metadata';
import { InputSection } from '@/sidebar/ui/InputSection';
import {
  droughtOptionsVariableState,
  droughtRcpParamState,
  droughtRiskVariableState,
  droughtShowOptionsState,
  droughtShowRiskState,
} from '@/state/drought/drought-parameters';

export const DroughtsControl: FC<{}> = () => {
  const [rcp, setRcp] = useRecoilState(droughtRcpParamState);
  const [showRisk, setShowRisk] = useRecoilState(droughtShowRiskState);
  const [showOptions, setShowOptions] = useRecoilState(droughtShowOptionsState);

  const [riskVariable, setRiskVariable] = useRecoilState(droughtRiskVariableState);
  const [optionsVariable, setOptionsVariable] = useRecoilState(droughtOptionsVariableState);

  return (
    <>
      <InputSection>
        <ParamDropdown<string>
          title="Climate Scenario (RCP)"
          value={rcp}
          onChange={setRcp}
          options={['2.6', '4.5', '8.5']}
        />
      </InputSection>
      <InputSection>
        <FormControlLabel
          label="Drought Risk"
          control={<Checkbox checked={showRisk} onChange={(e, checked) => setShowRisk(checked)} />}
        />
        <ParamDropdown<DroughtRiskVariableType>
          title="Drought Risk Variable"
          value={riskVariable}
          onChange={setRiskVariable}
          options={DROUGHT_RISK_VARIABLE_LABELS}
          disabled={!showRisk}
        />
        <Divider />
      </InputSection>
      <InputSection>
        <FormControlLabel
          label="Adaptation Options"
          control={<Checkbox checked={showOptions} onChange={(e, checked) => setShowOptions(checked)} />}
        />
        <ParamDropdown<DroughtOptionsVariableType>
          title="Adaptation Options Variable"
          value={optionsVariable}
          onChange={setOptionsVariable}
          options={DROUGHT_OPTIONS_VARIABLE_LABELS}
          disabled={!showOptions}
        />
      </InputSection>
    </>
  );
};
