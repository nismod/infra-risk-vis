import { Stack, Typography } from '@mui/material';
import { SidePanel } from 'details/SidePanel';
import { FC } from 'react';
import { FeatureAdaptationsTable } from './FeatureAdaptationsTable';

export const AdaptationsSidebar: FC<{}> = () => {
  return (
    <SidePanel height="80vh" px={2}>
      <Stack maxHeight="100%">
        <Typography variant="h6" gutterBottom>
          Adaptation Options
        </Typography>
        <FeatureAdaptationsTable />
      </Stack>
    </SidePanel>
  );
};
