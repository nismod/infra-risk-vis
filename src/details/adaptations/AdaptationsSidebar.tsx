import { Stack, Typography } from '@mui/material';
import { SidePanel } from 'details/SidePanel';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { FC } from 'react';
import { FeatureAdaptationsTable } from './FeatureAdaptationsTable';

export const AdaptationsSidebar: FC<{}> = () => {
  return (
    <SidePanel height="80vh" px={2}>
      <ErrorBoundary message="There was a problem displaying these details.">
        <Stack maxHeight="100%">
          <Typography variant="h6" gutterBottom>
            Adaptation Options
          </Typography>
          <FeatureAdaptationsTable />
        </Stack>
      </ErrorBoundary>
    </SidePanel>
  );
};
