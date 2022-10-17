import { FC } from 'react';
import { Box } from '@mui/material';

import { MapView } from '../map/MapView';
import { FeatureSidebar } from '../details/features/FeatureSidebar';
import { RegionDetails } from '../details/regions/RegionDetails';
import { SidebarContent } from 'sidebar/SidebarContent';
import { globalStyleVariables } from 'theme';
import { useSyncRecoilState } from 'lib/recoil/sync-state';
import { viewState, viewStateEffect } from 'state/view';
import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';
import { sectionStyleValueState, sectionVisibilityState } from 'state/sections';
import { selector, useRecoilValue } from 'recoil';
import { AdaptationsSidebar } from 'details/adaptations/AdaptationsSidebar';
import { SolutionsSidebar } from 'details/solutions/SolutionsSidebar';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';

interface MapViewProps {
  view: string;
}

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

const showAdaptationsTableState = selector<boolean>({
  key: 'showAdaptationsTable',
  get: ({ get }) => get(sectionVisibilityState('assets')) && get(sectionStyleValueState('assets')) === 'adaptation',
});

export const MapPage: FC<MapViewProps> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  const showAdaptationsTable = useRecoilValue(showAdaptationsTableState);
  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <SidebarLayout top={0} left={0} bottom={0} right={undefined} width={globalStyleVariables.controlSidebarWidth}>
        <ErrorBoundary message="There was a problem displaying the sidebar.">
          <SidebarContent />
        </ErrorBoundary>
      </SidebarLayout>
      <Box position="absolute" overflow="clip" top={globalStyleVariables.navbarHeight} left={0} right={0} bottom={0}>
        <ErrorBoundary message="There was a problem displaying the map." justifyErrorContent="center">
          <MapView />
        </ErrorBoundary>
      </Box>
      <SidebarLayout top={0} left={undefined} bottom={70} right={70} width={globalStyleVariables.detailSidebarWidth}>
        <Box mb={2}>
          <SolutionsSidebar />
        </Box>
        <Box mb={2}>{showAdaptationsTable ? <AdaptationsSidebar /> : <FeatureSidebar />}</Box>
        <Box mb={2}>
          <RegionDetails />
        </Box>
      </SidebarLayout>
    </ErrorBoundary>
  );
};
