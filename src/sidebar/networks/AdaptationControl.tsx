import { FormLabel } from '@mui/material';
import { CustomNumberSlider } from 'lib/controls/CustomSlider';
import { ParamDropdown } from 'lib/controls/ParamDropdown';
import { FC } from 'react';
import { useRecoilState } from 'recoil';
import { InputRow } from 'sidebar/ui/InputRow';
import { InputSection } from 'sidebar/ui/InputSection';
import { LayerStylePanel } from 'sidebar/ui/LayerStylePanel';
import { DataParam } from 'sidebar/ui/params/DataParam';
import { adaptationFieldState } from 'state/layers/networks';

function hazardLabel(val) {
  switch (val) {
    case 'flooding':
      return 'Flooding';
    case 'TC':
      return 'Cyclones';
    default:
      throw new Error('Unsupported hazard type: ' + val);
  }
}

function makeOptions(values, labelFn = (x) => x) {
  return values.map((val) => ({
    value: val,
    label: labelFn(val),
  }));
}

export const AdaptationControl: FC<{}> = () => {
  const [adaptationField, setAdaptationField] = useRecoilState(adaptationFieldState);
  return (
    <LayerStylePanel>
      <InputSection>
        <InputRow>
          <DataParam group="adaptation" id="sector">
            {({ value, onChange, options }) => (
              <ParamDropdown title="Sector" value={value} onChange={onChange} options={options} />
            )}
          </DataParam>
          <DataParam group="adaptation" id="subsector">
            {({ value, onChange, options }) => (
              <ParamDropdown title="Sub-sector" value={value} onChange={onChange} options={options} />
            )}
          </DataParam>
        </InputRow>
      </InputSection>
      <InputSection>
        <FormLabel>Adaptation for</FormLabel>
        <InputRow>
          <DataParam group="adaptation" id="hazard">
            {({ value, onChange, options }) => (
              <ParamDropdown
                title="Hazard"
                value={value}
                onChange={onChange}
                options={makeOptions(options, hazardLabel)}
              />
            )}
          </DataParam>
          <DataParam group="adaptation" id="rcp">
            {({ value, onChange, options }) => (
              <ParamDropdown title="RCP" value={value} onChange={onChange} options={options} />
            )}
          </DataParam>
        </InputRow>
      </InputSection>
      <InputSection>
        <DataParam group="adaptation" id="adaptation_name">
          {({ value, onChange, options }) => (
            <ParamDropdown title="Adaptation type" value={value} onChange={onChange} options={options} />
          )}
        </DataParam>
      </InputSection>
      <InputSection>
        <DataParam group="adaptation" id="adaptation_protection_level">
          {({ value, onChange, options }) =>
            options.length > 2 ? (
              <>
                <FormLabel>Protection level</FormLabel>
                <CustomNumberSlider title="Protection level" value={value} onChange={onChange} marks={options} />
              </>
            ) : (
              <ParamDropdown title="Protection level" value={value} onChange={onChange} options={options} />
            )
          }
        </DataParam>
      </InputSection>

      <ParamDropdown
        title="Displayed variable"
        value={adaptationField}
        onChange={setAdaptationField}
        options={[
          { value: 'avoided_ead_mean', label: 'Avoided Expected Annual Damages' },
          { value: 'avoided_eael_mean', label: 'Avoided Expected Annual Economic Losses' },
          { value: 'adaptation_cost', label: 'Adaptation Cost' },
        ]}
      />
    </LayerStylePanel>
  );
};
