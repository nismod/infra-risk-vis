import { FC } from 'react';
import { Box } from '@mui/material';

import { MapView } from './map/MapView';
import { SidebarContent } from 'sidebar/SidebarContent';
import { globalStyleVariables } from 'theme';

interface MapViewProps {
  view: string;
}

export const MapPage: FC<MapViewProps> = ({ view }) => {
  return (
    <>
      <Box
        position="absolute"
        top={globalStyleVariables.navbarHeight}
        left={0}
        bottom={0}
        width={globalStyleVariables.sidebarWidth}
        zIndex={1000}
        overflow="auto"
        boxSizing="border-box"
        my={1}
        sx={{ pointerEvents: 'none' }}
      >
        <Box px={1}>
          {/* <Toolbar /> */}
          {/* Prevents app bar from concealing content*/}
          <SidebarContent view={view} />
        </Box>
      </Box>
      <Box position="absolute" top={globalStyleVariables.navbarHeight} left={0} right={0} bottom={0}>
        <MapView view={view} />
      </Box>
    </>
  );
};
