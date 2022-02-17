import { ToggleButton, ToggleButtonGroup } from '@mui/material';
import { useRecoilState } from 'recoil';
import { viewModeState } from 'state/view-mode';

export const ViewModeToggle = () => {
  const [mode, setMode] = useRecoilState(viewModeState);

  return (
    <ToggleButtonGroup
      size="small"
      exclusive
      value={mode}
      onChange={(e, value) => setMode(value)}
      style={{ display: 'flex', flexDirection: 'row', justifyItems: 'stretch' }}
    >
      <ToggleButton style={{ width: '50%' }} value="input">
        Input Data
      </ToggleButton>
      <ToggleButton style={{ width: '50%' }} value="direct-damages">
        Direct Damages
      </ToggleButton>
    </ToggleButtonGroup>
  );
};
