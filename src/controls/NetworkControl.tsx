import React, { FC, useMemo } from 'react';
import { Box, Checkbox, FormControl, FormControlLabel, FormGroup, Typography } from '@material-ui/core';

import { LayerName, LAYERS } from '../config/layers';

interface NetworkControlProps {
  networkSelection: Record<string, boolean>;
  onNetworkSelection: (visUpdate: Record<string, boolean>) => void; //TODO change record key type to LayerName
}

export const NetworkControl: FC<NetworkControlProps> = ({ networkSelection, onNetworkSelection }) => {
  const networkLayers = useMemo(() => Object.keys(networkSelection), [networkSelection]);

  return (
    <Box mb={2}>
      <Typography variant="h6">Infrastructure Layers</Typography>
      <FormControl component="fieldset">
        {networkLayers.map((layerName) => {
          const { label, type, color } = LAYERS[layerName];
          const checked = networkSelection[layerName];
          return (
            <FormGroup row key={'toggleLayer' + layerName}>
              <FormControlLabel
                control={
                  <Checkbox
                    data-layer={layerName}
                    color="primary"
                    checked={checked}
                    value={layerName}
                    name={'toggleLayerCheckbox' + layerName}
                    onChange={(e) =>
                      onNetworkSelection({ ...networkSelection, [e.target.value as LayerName]: e.target.checked })
                    }
                  />
                }
                label={
                  <>
                    <span className={type === 'line' ? 'dot line' : 'dot'} style={{ backgroundColor: color }}></span>
                    {label}
                  </>
                }
              ></FormControlLabel>
            </FormGroup>
          );
        })}
      </FormControl>
    </Box>
  );
};
