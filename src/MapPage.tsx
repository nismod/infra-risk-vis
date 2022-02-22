import { FC } from 'react';
import { Box, Drawer, Toolbar } from '@mui/material';

import { FeatureSidebar } from './features/FeatureSidebar';
import { TooltipContent } from './map/tooltip/TooltipContent';
import { MapLegend } from './map/legend/MapLegend';
import { MapLayerSelection } from './map/layers/MapLayerSelection';
import { LegendContent } from './map/legend/LegendContent';

import { DataMapTooltip } from 'lib/data-map/DataMapTooltip';
import { MapSearch } from 'lib/map/place-search/MapSearch';

import { MapView } from './map/MapView';
import { SidebarContent } from 'sidebar/SidebarContent';

interface MapViewProps {
  view: string;
}

const sidebarWidth = 400;
const navBarHeight = 64;

export const MapPage: FC<MapViewProps> = ({ view }) => {
  return (
    <>
      <Box
        position="absolute"
        top={navBarHeight}
        left={0}
        bottom={0}
        width={sidebarWidth}
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
      <Box position="absolute" top={navBarHeight} left={0} right={0} bottom={0}>
        <MapView view={view}>
          <DataMapTooltip>
            <TooltipContent />
          </DataMapTooltip>
          <Box position="absolute" top={0} left={sidebarWidth + 10} ml={3} m={1} zIndex={1000}>
            <Box mb={1}>
              <MapSearch />
            </Box>
            <Box mb={1}>
              <MapLayerSelection />
            </Box>
          </Box>
          <Box position="absolute" bottom={0} left={sidebarWidth + 10} m={1} ml={3} zIndex={1000}>
            <MapLegend>
              <LegendContent />
            </MapLegend>
          </Box>
          <FeatureSidebar />
        </MapView>
      </Box>
    </>
  );
};
