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

  return (
    <Paper>
      <Box p={3}>
        <FeatureSidebarContent feature={feature} viewLayer={viewLayer} />
      </Box>
    </Paper>
  );
};
