import { TabContext, useTabContext } from '@mui/lab';
import { BottomNavigation, BottomNavigationAction, Box } from '@mui/material';
import React, { FC, useRef, useState } from 'react';
import { BottomSheet, BottomSheetRef } from 'react-spring-bottom-sheet';
import { useRecoilValue } from 'recoil';

import { MapView } from '@/map/MapView';
import { globalStyleVariables } from '@/theme';

import { mobileTabHasContentState } from './tab-has-content';
import { TabConfig, mobileTabsConfig } from './tabs-config';

import './bottom-sheet.css';

/**
 * Custom BottomNavigationAction that gets disabled if the corresponding tab doesn't have any content.
 *
 */
const TabNavigationAction: FC<{
  value: string;

  label: TabConfig['label'];
  IconComponent: TabConfig['IconComponent'];

  selected?: boolean;
  showLabel?: boolean;
  onChange?: any;
}> = ({ value, label, IconComponent, selected, showLabel, onChange }) => {
  const hasContent = useRecoilValue(mobileTabHasContentState(value));
  const disabled = !hasContent;

  // cloneElement is needed here because MUI BottomNavigation uses Children.map and cloneElement
  // and expects the children to be BottomNavigationActions, so this component needs to adapt to that
  // see @mui/material/BottomNavigation/BottomNavigation.js#L64
  return React.cloneElement(
    <BottomNavigationAction
      value={value}
      label={label}
      icon={<IconComponent fontSize="small" />}
      disabled={disabled}
      sx={
        selected
          ? undefined
          : {
              '&.Mui-disabled': {
                color: '#b2b2b2',
              },
            }
      }
      showLabel={showLabel}
      onChange={onChange}
    />,
    {
      selected,
    },
  );
};

/**
 * TabPanel that doesn't unmount inactive tabs, only uses display:none
 */
const NonUnmountingTabPanel = ({ value, children }) => {
  const context = useTabContext();

  return <Box display={value === context?.value ? 'block' : 'none'}>{children}</Box>;
};

const MobileTabPanel: FC<{ tabConfig: TabConfig }> = ({ tabConfig: { id, ContentComponent } }) => (
  <NonUnmountingTabPanel value={id}>
    <Box m={2}>
      <ContentComponent />
    </Box>
  </NonUnmountingTabPanel>
);

export const MapPageMobileLayout = () => {
  const [bottomTabId, setBottomTabId] = useState('layers');
  const sheetRef = useRef<BottomSheetRef>();

  /*
    const handleChange = useCallback(() => {
      if (sheetRef.current.height <= 100) {
        sheetRef.current.snapTo(({ snapPoints }) => snapPoints[1]);
      }
    }, []);
    */

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
          <Box sx={{ cursor: 'default' }}>
            <BottomNavigation
              showLabels
              value={bottomTabId}
              onChange={(e, value) => setBottomTabId(value)}
            >
              {mobileTabsConfig.map(({ id, label, IconComponent }) => {
                return (
                  <TabNavigationAction
                    key={id}
                    label={label}
                    IconComponent={IconComponent}
                    value={id}
                  />
                );
              })}
            </BottomNavigation>
          </Box>
        }
        header={<Box height={10} />}
      >
        <TabContext value={bottomTabId}>
          {mobileTabsConfig.map((tabConfig) => (
            <MobileTabPanel key={tabConfig.id} tabConfig={tabConfig} />
          ))}
        </TabContext>
      </BottomSheet>
    </>
  );
};
