import { useRecoilValue } from 'recoil';

import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { Box } from '@mui/material';
import { RegionDetailsContent } from './RegionDetailsContent';
import { SidePanel } from 'details/SidePanel';
import { DeselectButton } from 'details/DeselectButton';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';

export const RegionDetails = () => {
  const selectedRegion = useRecoilValue(selectionState('regions'));

  if (!selectedRegion) return null;

  return (
    <SidePanel position="relative">
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box position="absolute" top={0} right={0} p={2}>
          <DeselectButton interactionGroup="regions" title="Deselect region" />
        </Box>
        <RegionDetailsContent selectedRegion={selectedRegion} />
      </ErrorBoundary>
    </SidePanel>
  );
};
