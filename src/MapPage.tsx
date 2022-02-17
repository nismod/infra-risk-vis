import { FC } from 'react';
import { Box, Drawer, Toolbar } from '@mui/material';

import { DataMap } from './map/DataMap';
import { ViewName } from './config/views';
import { SidebarContent } from 'sidebar/SidebarContent';

interface MapViewProps {
  view: ViewName;
}

const sidebarWidth = 360;

export const MapPage: FC<MapViewProps> = ({ view }) => {
  return (
    <>
      <Drawer variant="permanent">
        <Box p={4} pt={2} width={sidebarWidth} boxSizing="border-box">
          <Toolbar /> {/* Prevents app bar from concealing content*/}
          <SidebarContent view={view} />
        </Box>
      </Drawer>
      <Box position="absolute" top={64} left={sidebarWidth} right={0} bottom={0}>
        <DataMap view={view} />
      </Box>
    </>
  );
};
