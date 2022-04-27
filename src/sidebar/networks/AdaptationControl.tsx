import { ParamDropdown } from 'lib/controls/ParamDropdown';
import { FC } from 'react';
import { useRecoilState } from 'recoil';
import { LayerStylePanel } from 'sidebar/ui/LayerStylePanel';
import { adaptationFieldState } from 'state/layers/networks';

export const AdaptationControl: FC<{}> = () => {
  const [adaptationField, setAdaptationField] = useRecoilState(adaptationFieldState);
  return (
    <LayerStylePanel>
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
