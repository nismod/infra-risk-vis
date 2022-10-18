import { Box } from '@mui/system';
import React, { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { selectionState } from '@/lib/data-map/interactions/interaction-state';
import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { DeselectButton } from '@/details/DeselectButton';
import { SidePanel } from '@/details/SidePanel';

import { FeatureSidebarContent } from './FeatureSidebarContent';

export const FeatureSidebar: FC<{}> = () => {
  const featureSelection = useRecoilValue(selectionState('assets'));

  if (!featureSelection) return null;

  const {
    target: { feature },
    viewLayer,
  } = featureSelection;

  return (
    <SidePanel position="relative">
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box position="absolute" top={15} right={15} zIndex={1000}>
          <DeselectButton interactionGroup="assets" title="Deselect asset" />
        </Box>
        <FeatureSidebarContent feature={feature} assetType={viewLayer.id} />
      </ErrorBoundary>
    </SidePanel>
  );
};
