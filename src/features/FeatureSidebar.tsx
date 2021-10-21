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
    <Box position="absolute" top={0} right={0} width={300} m={2}>
      <Paper>
        <Box p={3} overflow="auto" height="90vh">
          <FeatureSidebarContent f={f} layer={logicalLayerConfig} />
        </Box>
      </Paper>
    </Box>
  );
};
