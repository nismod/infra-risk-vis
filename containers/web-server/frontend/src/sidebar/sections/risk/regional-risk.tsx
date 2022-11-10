import { useRecoilState } from 'recoil';

import { ParamDropdown } from '@/lib/controls/ParamDropdown';

import { REGIONAL_EXPOSURE_VARIABLE_LABELS, RegionalExposureVariableType } from '@/config/regional-risk/metadata';
import { InputSection } from '@/sidebar/ui/InputSection';
import { regionalExposureVariableState } from '@/state/data-selection/regional-risk';


export const RegionalRiskSection = () => {
    const [rexpVariable, setRexpVariable] = useRecoilState(regionalExposureVariableState);

    return (
      <>
        <InputSection>
          <ParamDropdown<RegionalExposureVariableType>
            title="Population exposed to:"
            value={rexpVariable}
            onChange={setRexpVariable}
            options={REGIONAL_EXPOSURE_VARIABLE_LABELS}
          />
        </InputSection>
      </>
    )
}
