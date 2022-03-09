import { FC } from 'react';
import { Box } from '@mui/material';

import { MapView } from './map/MapView';
import { SidebarContent } from 'sidebar/SidebarContent';
import { globalStyleVariables } from 'theme';
import { useSyncRecoilState } from 'lib/recoil/sync-state';
import { viewState } from 'state/view';

interface MapViewProps {
  view: string;
}

const SidebarLayout = ({ children }) => (
  <Box
    position="absolute"
    top={globalStyleVariables.navbarHeight}
    left={0}
    bottom={0}
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
      <SidebarLayout>
        <SidebarContent />
      </SidebarLayout>
      <Box position="absolute" overflow="clip" top={globalStyleVariables.navbarHeight} left={0} right={0} bottom={0}>
        <MapView />
      </Box>
    </>
  );
};
