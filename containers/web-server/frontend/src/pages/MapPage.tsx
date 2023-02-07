import { Layers, Palette, TableRows } from '@mui/icons-material';
import { BottomNavigation, BottomNavigationAction, Box } from '@mui/material';
import { FC, useRef, useState } from 'react';
import { BottomSheet, BottomSheetRef } from 'react-spring-bottom-sheet';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';
import { useSyncRecoilState } from '@/lib/recoil/sync-state';

import { InitData } from '@/InitData';
import { DetailsContent } from '@/details/DetailsContent';
import { MapView } from '@/map/MapView';
import { MapLegend } from '@/map/legend/MapLegend';
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

const MapPageMobileLayout = () => {
  const [bottomTab, setBottomTab] = useState(0);
  const sheetRef = useRef<BottomSheetRef>();

  /*
  const handleChange = useCallback(() => {
    if (sheetRef.current.height <= 100) {
      sheetRef.current.snapTo(({ snapPoints }) => snapPoints[1]);
    }
  }, []);
  */

  const tabs = [
    {
      label: 'Layers',
      icon: <Layers fontSize="small" />,
      content: <LayersSidebar key="layers" />,
    },
    {
      label: 'Legend',
      icon: <Palette fontSize="small" />,
      content: <MapLegend key="legend" />,
    },
    {
      label: 'Selection',
      icon: <TableRows fontSize="small" />,
      content: <DetailsContent key="details" />,
    },
  ];

  return (
    <>
      <Box position="absolute" overflow="clip" top={0} left={0} right={0} bottom={0}>
        <MapView />
      </Box>
      <BottomSheet
        ref={sheetRef}
        open={true}
        blocking={false}
        expandOnContentDrag={true}
        skipInitialTransition={true}
        snapPoints={({ footerHeight, maxHeight }) => [
          footerHeight + 38, // magic number to allow the puller/header to be visible at the smallest snap point
          maxHeight * 0.35,
          maxHeight * 0.6,
          maxHeight - globalStyleVariables.navbarHeight - 4,
        ]}
        footer={
          <BottomNavigation
            showLabels
            // sx={{ padding: 0 }}
            value={bottomTab}
            onChange={(e, value) => setBottomTab(value)}
          >
            {tabs.map(({ label, icon }) => (
              <BottomNavigationAction key={label} label={label} icon={icon} />
            ))}
          </BottomNavigation>
        }
        header={<Box height={10} />}
      >
        <Box m={2}>{tabs[bottomTab].content}</Box>
      </BottomSheet>
    </>
  );
};

const MapPageLayout = () => {
  const isMobile = useIsMobile();

  return isMobile ? <MapPageMobileLayout /> : <MapPageDesktopLayout />;
};

export const MapPage: FC<{ view: string }> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <InitData />
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <MapPageLayout />
    </ErrorBoundary>
  );
};
