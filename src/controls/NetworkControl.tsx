import React, { FC, useMemo } from 'react';
import { Box, Checkbox, FormControl, FormControlLabel, FormGroup, Typography } from '@material-ui/core';

import { LayerName, LAYERS } from '../config/layers';

interface NetworkControlProps {
  networkSelection: Record<string, boolean>;
  onNetworkSelection: (visUpdate: Record<string, boolean>) => void; //TODO change record key type to LayerName
}

const LayerLabel = ({ layerConfig: { label, type, color } }) => {
  return (
    <Typography>
      <Box className={type === 'line' ? 'dot line' : 'dot'} style={{ backgroundColor: color }}></Box>
      {label}
    </Typography>
  );
};
export const NetworkControl: FC<NetworkControlProps> = ({ networkSelection, onNetworkSelection }) => {
  const networkLayers = useMemo(() => Object.keys(networkSelection), [networkSelection]);

  return (
    <Box mb={2}>
      <Typography variant="h6">Infrastructure Assets</Typography>
      <FormControl component="fieldset">
        {networkLayers.map((layerName) => (
          <FormControlLabel
            key={layerName}
            style={{ marginBottom: -12 }}
            control={
              <Checkbox
                data-layer={layerName}
                color="default"
                checked={networkSelection[layerName]}
                value={layerName}
                name={'toggleLayerCheckbox' + layerName}
                onChange={(e) =>
                  onNetworkSelection({ ...networkSelection, [e.target.value as LayerName]: e.target.checked })
                }
                // style={{ padding: 4 }}
              />
            }
            label={<LayerLabel layerConfig={LAYERS[layerName]} />}
          ></FormControlLabel>
        ))}
      </FormControl>
    </Box>
  );
};
