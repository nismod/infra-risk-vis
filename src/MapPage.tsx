import { FC } from 'react';
import { Box } from '@mui/material';

import { MapView } from './map/MapView';
import { FeatureSidebar } from './details/features/FeatureSidebar';
import { RegionDetails } from './details/regions/RegionDetails';
import { SidebarContent } from 'sidebar/SidebarContent';
import { globalStyleVariables } from 'theme';
import { useSyncRecoilState } from 'lib/recoil/sync-state';
import { viewState, viewStateEffect } from 'state/view';
import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';

interface MapViewProps {
  view: string;
}

const SidebarLayout = ({ top, bottom, left, right, children }) => (
  <Box
    position="absolute"
    top={globalStyleVariables.navbarHeight + top}
    bottom={bottom}
    left={left}
    right={right}
    width={globalStyleVariables.sidebarWidth - 20}
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
      <Box width={globalStyleVariables.sidebarWidth - 20}>
        <Box pl={1} pt={1} sx={{ pointerEvents: 'auto' }}>
          {children}
        </Box>
      </Box>
    </Box>
  </Box>
);

export const MapPage: FC<MapViewProps> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  return (
    <>
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <SidebarLayout top={0} left={0} bottom={0} right={undefined}>
        <SidebarContent />
      </SidebarLayout>
      <Box position="absolute" overflow="clip" top={globalStyleVariables.navbarHeight} left={0} right={0} bottom={0}>
        <MapView />
      </Box>
      <SidebarLayout top={0} left={undefined} bottom={70} right={70}>
        <Box mb={2}>
          <FeatureSidebar />
        </Box>
        <Box mb={2}>
          <RegionDetails />
        </Box>
      </SidebarLayout>
    </>
  );
};
