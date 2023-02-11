import { Layers, Palette, TableRows } from '@mui/icons-material';
import { BottomNavigation, BottomNavigationAction, Box } from '@mui/material';
import { useRef, useState } from 'react';
import { BottomSheet, BottomSheetRef } from 'react-spring-bottom-sheet';

import { DetailsContent } from '@/details/DetailsContent';
import { MapView } from '@/map/MapView';
import { MapLegend } from '@/map/legend/MapLegend';
import { LayersSidebar } from '@/sidebar/LayersSidebar';
import { globalStyleVariables } from '@/theme';

export const MapPageMobileLayout = () => {
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
