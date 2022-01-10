import React, { FC } from 'react';
import { Box, Paper } from '@material-ui/core';

import { VectorHover } from '../map/DataMap';
import { LAYERS } from '../config/layers';
import { FeatureSidebarContent } from './FeatureSidebarContent';

interface FeatureSidebarProps {
  featureSelection: VectorHover;
}

export const FeatureSidebar: FC<FeatureSidebarProps> = ({ featureSelection }) => {
  const { feature, logicalLayer } = featureSelection;
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
