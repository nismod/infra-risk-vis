import { ToggleSection, ToggleSectionGroup } from 'lib/controls/accordion-toggle/ToggleSection';

import { hazardSelectionState } from '../../state/hazards/hazard-selection';
import { InputRow } from 'sidebar/ui/InputRow';
import { InputSection } from 'sidebar/ui/InputSection';
import { ReturnPeriodControl } from 'sidebar/ui/params/ReturnPeriodControl';
import { EpochControl } from 'sidebar/ui/params/EpochControl';
import { RCPControl } from 'sidebar/ui/params/RCPControl';

export const HazardsControl = () => {
  const forceSingle = false;

  return (
    <>
      <ToggleSectionGroup toggleState={hazardSelectionState}>
        <ToggleSection id="fluvial" label="River Flooding" forceSingle={forceSingle}>
          <ReturnPeriodControl group="fluvial" param="returnPeriod" />
        </ToggleSection>

        <ToggleSection id="surface" label="Surface Flooding" forceSingle={forceSingle}>
          <ReturnPeriodControl group="surface" param="returnPeriod" />
        </ToggleSection>

        <ToggleSection id="coastal" label="Coastal Flooding" forceSingle={forceSingle}>
          <InputSection>
            <ReturnPeriodControl group="coastal" param="returnPeriod" />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="coastal" />
              <RCPControl group="coastal" />
            </InputRow>
          </InputSection>
        </ToggleSection>

        <ToggleSection id="cyclone" label="Cyclones" forceSingle={forceSingle}>
          <InputSection>
            <ReturnPeriodControl
              group="cyclone"
              valueLabelDisplay="auto"
              showMarkLabelsFor={[10, 50, 100, 500, 1000, 5000, 10000]}
            />
          </InputSection>
          <InputSection>
            <InputRow>
              <EpochControl group="cyclone" />
              <RCPControl group="cyclone" />
            </InputRow>
          </InputSection>
        </ToggleSection>
      </ToggleSectionGroup>
    </>
  );
};
