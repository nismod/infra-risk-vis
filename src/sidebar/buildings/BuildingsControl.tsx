import { Checkbox, FormControlLabel, FormGroup } from '@mui/material';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { useRecoilState } from 'recoil';
import { LayerLabel } from 'sidebar/ui/LayerLabel';
import { buildingSelectionState } from 'state/buildings';

export const BuildingsControl = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(buildingSelectionState);
  return (
    <FormGroup>
      {Object.entries(checkboxState).map(([key, value]) => (
        <FormControlLabel
          key={key}
          control={
            <Checkbox
              checked={value}
              onChange={(e, checked) => setCheckboxState({ ...checkboxState, [key]: checked })}
            />
          }
          label={<LayerLabel {...NETWORKS_METADATA[key]} />}
        />
      ))}
    </FormGroup>
  );
};
