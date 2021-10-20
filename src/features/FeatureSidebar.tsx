import React, { FC } from 'react';
import { Box, List, ListItem, ListItemText, Paper, Typography } from '@material-ui/core';

import { titleCase } from '../helpers';
import { VectorHover } from '../map/DataMap';
import { LAYERS } from '../config/layers';

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
          <Typography variant="h6">Selected Asset</Typography>
          <Typography variant="body1" style={{ color: logicalLayerConfig.color }}>
            {logicalLayerConfig.label}
          </Typography>
          <pre style={{ display: 'none' }}>
            <code>{JSON.stringify(f, null, 2)}</code>
          </pre>
          <List>
            {Object.entries(f).map(([key, value]) => (
              <ListItem key={key}>
                <ListItemText
                  primary={titleCase(key.replace(/_/g, ' '))}
                  primaryTypographyProps={{ variant: 'caption' }}
                  secondary={value === ' ' ? '-' : value || '-'}
                />
              </ListItem>
            ))}
          </List>
        </Box>
      </Paper>
    </Box>
  );
};
