import { ToggleButtonGroup, ToggleButton } from '@mui/material';
import { useRecoilState } from 'recoil';
import { RegionLevel } from 'config/regions/metadata';
import { regionLevelState } from 'state/regions';
import { useCallback } from 'react';

export const RegionLevelSelection = () => {
  const [regionLevel, setRegionLevel] = useRecoilState(regionLevelState);

  const handleChange = useCallback(
    (e, value: string) => {
      if (value != null) {
        setRegionLevel(value as RegionLevel);
      }
    },
    [setRegionLevel],
  );

  return (
    <ToggleButtonGroup exclusive value={regionLevel} onChange={handleChange} fullWidth>
      <ToggleButton value="level0" sx={{ textTransform: 'none' }}>
        Countries
      </ToggleButton>
      <ToggleButton value="level1" sx={{ textTransform: 'none' }}>
        Regions
      </ToggleButton>
    </ToggleButtonGroup>
  );
};

export const RegionsControl = () => {
  return <RegionLevelSelection />;
};
