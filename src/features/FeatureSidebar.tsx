import React, { FC } from 'react';
import { Box, Paper } from '@mui/material';

import { FeatureSidebarContent } from './FeatureSidebarContent';
import { useRecoilValue } from 'recoil';
import { selectionState } from 'lib/data-map/interactions/interaction-state';

export const FeatureSidebar: FC<{}> = () => {
  const featureSelection = useRecoilValue(selectionState('assets'));

  if (!featureSelection) return null;

  const {
    target: { feature },
    viewLayer,
  } = featureSelection;
  const f = feature.properties;

  return (
    <Box position="absolute" top={10} right={50} width={400}>
      <Paper>
        <Box p={3} maxHeight="calc(100vh - 125px)" style={{ overflowY: 'scroll' }}>
          <FeatureSidebarContent f={f} viewLayer={viewLayer} />
        </Box>
      </Paper>
    </Box>
  );
};
