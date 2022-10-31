import { Box } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { selectionState } from '@/lib/data-map/interactions/interaction-state';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { DeselectButton } from '@/details/DeselectButton';
import { SidePanel } from '@/details/SidePanel';

import { SolutionsSidebarContent } from './SolutionsSidebarContent';

export const SolutionsSidebar: FC<{}> = () => {
  const featureSelection = useRecoilValue(selectionState('solutions'));

  if (!featureSelection) return null;

  const {
    target: { feature },
    viewLayer,
  } = featureSelection;

  return (
    <SidePanel position="relative">
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box position="absolute" top={0} right={0} p={2}>
          <DeselectButton interactionGroup="solutions" title="Deselect" />
        </Box>
        <SolutionsSidebarContent feature={feature} solutionType={viewLayer.id} />
      </ErrorBoundary>
    </SidePanel>
  );
};
