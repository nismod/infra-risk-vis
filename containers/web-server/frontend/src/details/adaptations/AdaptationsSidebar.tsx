import { Typography } from '@mui/material';
import { Box } from '@mui/system';
import { SidePanel } from 'details/SidePanel';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { FC } from 'react';
import { FeatureAdaptationsTable } from './FeatureAdaptationsTable';

export const AdaptationsSidebar: FC<{}> = () => {
  return (
    <SidePanel height="80vh" pb={1} px={0} pt={0}>
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box px={2} pt={2}>
          <Typography variant="h6" gutterBottom>
            Adaptation Options
          </Typography>
        </Box>
        <Box height="73vh" position="relative">
          <FeatureAdaptationsTable />
        </Box>
      </ErrorBoundary>
    </SidePanel>
  );
};
