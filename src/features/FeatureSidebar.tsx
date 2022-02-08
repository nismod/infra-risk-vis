import React, { FC } from 'react';
import { Box, Paper } from '@mui/material';

import { LAYERS } from '../config/layers';
import { FeatureSidebarContent } from './FeatureSidebarContent';
import { InteractionTarget } from 'lib/map/interactions/use-interactions';

interface FeatureSidebarProps {
  featureSelection: InteractionTarget<any>;
}

export const FeatureSidebar: FC<FeatureSidebarProps> = ({ featureSelection }) => {
  const { target: feature, logicalLayer } = featureSelection;
  const f = feature.properties;
  const logicalLayerConfig = LAYERS[logicalLayer];

  return (
    <Box position="absolute" top={10} right={50} width={400}>
      <Paper>
        <Box p={3} maxHeight="calc(100vh - 125px)" style={{ overflowY: 'scroll' }}>
          <FeatureSidebarContent f={f} layer={logicalLayerConfig} />
        </Box>
      </Paper>
    </Box>
  );
};
