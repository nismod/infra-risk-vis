import { Typography } from '@mui/material';
import { Box } from '@mui/system';
import { HelpNote } from 'assessment/HelpNote';
import { getAssetDataFormats } from 'config/assets/data-formats';
import { SidePanel } from 'details/SidePanel';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import _ from 'lodash';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';
import { adaptationFieldSpecState, adaptationLayerSpecState } from 'state/layers/networks';
import { FeatureAdaptationsTable } from './FeatureAdaptationsTable';

export const AdaptationsSidebar: FC<{}> = () => {
  const { sector, subsector, assetType } = useRecoilValue(adaptationLayerSpecState);
  const fieldSpec = useRecoilValue(adaptationFieldSpecState);
  const { getDataLabel } = getAssetDataFormats(fieldSpec);
  return (
    <SidePanel height="80vh" pb={1} px={0} pt={0}>
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box px={2} pt={2}>
          <Typography variant="h6" gutterBottom>
            Adaptation Options
          </Typography>
        </Box>
        <HelpNote sx={{ mx: 1 }}>
          This table shows all {_.startCase(subsector)} ({_.startCase(assetType)}) assets, sorted by{' '}
          {getDataLabel(fieldSpec)} in descending order.
        </HelpNote>
        <Box height="67vh" position="relative">
          <FeatureAdaptationsTable />
        </Box>
      </ErrorBoundary>
    </SidePanel>
  );
};
