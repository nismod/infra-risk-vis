import { ToggleSection, ToggleSectionGroup } from 'lib/controls/accordion-toggle/ToggleSection';

import { hazardSelectionState } from '../../state/hazards/hazard-selection';
import { InputRow } from 'sidebar/ui/InputRow';
import { InputSection } from 'sidebar/ui/InputSection';
import { ReturnPeriodControl } from 'sidebar/ui/params/ReturnPeriodControl';
import { EpochControl } from 'sidebar/ui/params/EpochControl';
import { RCPControl } from 'sidebar/ui/params/RCPControl';
import { useRecoilValue } from 'recoil';
import { showDamagesState } from 'state/damage-mapping/damage-map';
import { Alert, Box } from '@mui/material';

export const HazardsControl = () => {
  const showDirectDamages = useRecoilValue(showDamagesState);
  const disabled = showDirectDamages;

  return (
    <>
      {showDirectDamages ? (
        <Box my={1}>
          <Alert severity="info">Hazards are currently following the Infrastructure &gt; Damages selection</Alert>
        </Box>
      ) : null}
      <ToggleSectionGroup toggleState={hazardSelectionState}>
        <ToggleSection id="fluvial" label="River Flooding" disabled={disabled}>
          <InputSection>
            <ReturnPeriodControl group="fluvial" param="returnPeriod" disabled={disabled} />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="fluvial" disabled={disabled} />
              <RCPControl group="fluvial" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="coastal" label="Coastal Flooding" disabled={disabled}>
          <InputSection>
            <ReturnPeriodControl group="coastal" param="returnPeriod" disabled={disabled} />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="coastal" disabled={disabled} />
              <RCPControl group="coastal" disabled={disabled} />
            </InputRow>
          </InputSection>
        </ToggleSection>

      </ToggleSectionGroup>
    </>
  );
};
