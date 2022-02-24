import { FormControlLabel, Radio, RadioGroup, Box, FormControl, FormLabel, Select, MenuItem } from '@mui/material';
import { useRecoilState } from 'recoil';
import { RegionLevel } from 'config/regions/metadata';
import { regionDataState, RegionDataType, regionLevelState } from 'state/regions';

export const RegionLevelSelection = () => {
  const [boundaryLevel, setBoundaryLevel] = useRecoilState(regionLevelState);
  return (
    <FormControl component="fieldset">
      <FormLabel>Administrative division</FormLabel>
      <RadioGroup value={boundaryLevel} onChange={(e, value) => setBoundaryLevel(value as RegionLevel)}>
        <FormControlLabel value="parish" label="Parishes" control={<Radio unselectable="off" />} />
        <FormControlLabel value="enumeration" label="Enumeration Districts" control={<Radio unselectable="off" />} />
      </RadioGroup>
    </FormControl>
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
