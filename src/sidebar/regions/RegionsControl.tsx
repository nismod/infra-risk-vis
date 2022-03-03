import { Box, FormControl, FormLabel, Select, MenuItem, ToggleButtonGroup, ToggleButton } from '@mui/material';
import { useRecoilState } from 'recoil';
import { RegionLevel } from 'config/regions/metadata';
import { regionDataState, RegionDataType, regionLevelState } from 'state/regions';
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
      <ToggleButton value="parish" sx={{ textTransform: 'none' }}>
        Parishes
      </ToggleButton>
      <ToggleButton value="enumeration" sx={{ textTransform: 'none' }}>
        Enumeration Districts
      </ToggleButton>
    </ToggleButtonGroup>
  );
};

const RegionsDataSelect = () => {
  const [regionData, setRegionData] = useRecoilState(regionDataState);
  return (
    <FormControl fullWidth>
      <FormLabel>Mapped data</FormLabel>
      <Select variant="standard" value={regionData} onChange={(e) => setRegionData(e.target.value as RegionDataType)}>
        <MenuItem value="boundaries">None (boundaries only)</MenuItem>
        <MenuItem value="population">Population</MenuItem>
      </Select>
    </FormControl>
  );
};

export const RegionsControl = () => {
  return (
    <>
      <Box mb={2}>
        <RegionLevelSelection />
      </Box>
      <Box>
        <RegionsDataSelect />
      </Box>
    </>
  );
};
