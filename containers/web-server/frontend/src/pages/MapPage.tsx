import { Box } from '@mui/material';
import { FC } from 'react';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';
import { useSyncRecoilState } from '@/lib/recoil/sync-state';

import { InitData } from '@/InitData';
import { DetailsContent } from '@/details/DetailsContent';
import { MapView } from '@/map/MapView';
import { LayersSidebar } from '@/sidebar/LayersSidebar';
import { viewState, viewStateEffect } from '@/state/view';
import { globalStyleVariables } from '@/theme';
import { useIsMobile } from '@/use-is-mobile';

const SidebarLayout = ({ top, bottom, left, right, width, children }) => (
  <Box
    position="absolute"
    top={top}
    bottom={bottom}
    left={left}
    right={right}
    width={width - 20}
    zIndex={1000}
    sx={{ pointerEvents: 'none' }}
  >
    <Box
      overflow="auto"
      maxHeight="100%"
      sx={{ pointerEvents: 'auto' }}
      position="absolute"
      left={0}
      right={-20}
      top={0}
    >
      <Box width={width - 20}>
        <Box pl={1} pt={1} sx={{ pointerEvents: 'auto' }}>
          {children}
        </Box>
      </Box>
    </Box>
  </Box>
);

const MapPageDesktopLayout = () => {
  return (
    <>
      <SidebarLayout
        top={0}
        left={0}
        bottom={0}
        right={undefined}
        width={globalStyleVariables.controlSidebarWidth}
      >
        <LayersSidebar />
      </SidebarLayout>
      <Box position="absolute" overflow="clip" top={0} left={0} right={0} bottom={0}>
        <MapView />
      </Box>
      <SidebarLayout
        top={0}
        left={undefined}
        bottom={70}
        right={70}
        width={globalStyleVariables.detailSidebarWidth}
      >
        <DetailsContent />
      </SidebarLayout>
    </>
  );
};

export const MapPage: FC<{ view: string }> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <InitData />
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <MapPageDesktopLayout />
    </ErrorBoundary>
  );
};
