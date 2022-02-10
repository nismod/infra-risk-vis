import { FC, useMemo } from 'react';
import { Box, MenuItem, Select, Typography } from '@mui/material';

import { titleCase, unique } from 'lib/helpers';
import { useSelect } from 'lib/hooks/use-select';

import { EADChart } from './charts/EADChart';

interface EADDataPoint {
  hazardType: string;
}
export const EADChartSection: FC<{ eadData: EADDataPoint[] }> = ({ eadData }) => {
  const hazards = useMemo(() => unique(eadData.map((d) => d.hazardType)), [eadData]);
  const [selectedHazard, setSelectedHazard] = useSelect(hazards);

  const selectedData = useMemo(
    () => (selectedHazard ? eadData.filter((x) => x.hazardType === selectedHazard) : null),
    [selectedHazard, eadData],
  );
  return (
    <Box py={2}>
      <Box sx={{ justifyContent: 'space-between' }}>
        <Typography variant="subtitle2">Expected Annual Damages (EAD)</Typography>
        {hazards.length !== 0 && (
          <Select
            variant="standard"
            value={selectedHazard ?? ''}
            onChange={(e) => setSelectedHazard(e.target.value as string)}
          >
            {hazards.map((h) => (
              <MenuItem key={h} value={h}>
                {titleCase(h)}
              </MenuItem>
            ))}
          </Select>
        )}
      </Box>
      {selectedData ? (
        <Box mt={1}>
          <EADChart
            data={{
              table: selectedData,
            }}
            actions={false}
            padding={0}
            width={260} // this is currently picked to fit the chart to the sidebar width
            height={150}
            renderer="svg"
          />
        </Box>
      ) : (
        <Typography variant="body2" color="textSecondary">
          No exposure direct damages estimated.
        </Typography>
      )}
    </Box>
  );
};
