import { FC, useMemo } from 'react';
import { Box, Typography } from '@material-ui/core';

import { titleCase } from '../helpers';

import { EADChart } from './charts/EADChart';

export const EADChartSection: FC<{ eadData: any }> = ({ eadData }) => {
  // TODO enable hazard selection
  const selectedHazard = eadData.map((x) => x.hazardType)[0];
  const selectedData = useMemo(
    () => (selectedHazard ? eadData.filter((x) => x.hazardType === selectedHazard) : null),
    [selectedHazard, eadData],
  );
  return (
    <Box py={2}>
      <Typography variant="subtitle2">
        Expected Annual Damages (EAD) {selectedData && <>&ndash; {titleCase(selectedHazard)}</>}
      </Typography>
      {selectedData ? (
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
      ) : (
        <Typography variant="body2" color="textSecondary">
          No exposure direct damages estimated.
        </Typography>
      )}
    </Box>
  );
};
