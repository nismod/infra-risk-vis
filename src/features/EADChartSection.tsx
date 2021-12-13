import { FC } from 'react';
import { Box, Typography } from '@material-ui/core';

import { titleCase } from '../helpers';

import { EADChart } from './charts/EADChart';

export const EADChartSection: FC<{ eadData: any }> = ({ eadData }) => {
  // TODO enable hazard selection
  const selectedHazard = eadData.map((x) => x.hazardType)[0];
  const selectedData = selectedHazard ? eadData.filter((x) => x.hazardType === selectedHazard) : null;
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
          width={260}
        />
      ) : (
        <Typography variant="body2" color="textSecondary">
          No exposure direct damages estimated.
        </Typography>
      )}
    </Box>
  );
};
