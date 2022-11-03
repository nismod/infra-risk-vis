import { Box, Paper } from '@mui/material';
import { FC } from 'react';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';
import { useSyncRecoilState } from '@/lib/recoil/sync-state';

import { DetailsContent } from '@/details/DetailsContent';
import { MapView } from '@/map/MapView';
import { SidebarContent } from '@/sidebar/SidebarContent';
import { viewState, viewStateEffect } from '@/state/view';
import { globalStyleVariables } from '@/theme';

const SidebarLayout = ({ top, bottom, left, right, width, children }) => (
  <Box
    position="absolute"
    top={globalStyleVariables.navbarHeight + top}
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

interface MapPageProps {
  view: string;
}

export const MapPage: FC<MapPageProps> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <SidebarLayout top={0} left={0} bottom={0} right={undefined} width={globalStyleVariables.controlSidebarWidth}>
        <Paper elevation={0}>
          <ErrorBoundary message="There was a problem displaying the sidebar.">
            <SidebarContent />
          </ErrorBoundary>
        </Paper>
      </SidebarLayout>
      <Box position="absolute" overflow="clip" top={globalStyleVariables.navbarHeight} left={0} right={0} bottom={0}>
        <ErrorBoundary message="There was a problem displaying the map." justifyErrorContent="center">
          <MapView />
        </ErrorBoundary>
      </Box>
      <SidebarLayout top={0} left={undefined} bottom={70} right={70} width={globalStyleVariables.detailSidebarWidth}>
        <DetailsContent />
      </SidebarLayout>
    </ErrorBoundary>
  );
};
