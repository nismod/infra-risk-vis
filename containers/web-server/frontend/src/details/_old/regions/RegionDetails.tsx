import { Box } from '@mui/material';
import { useRecoilValue } from 'recoil';

import { selectionState } from '@/lib/data-map/interactions/interaction-state';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { DeselectButton } from '@/details/ui/DeselectButton';
import { SidePanel } from '@/details/ui/SidePanel';

import { RegionDetailsContent } from './RegionDetailsContent';

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
