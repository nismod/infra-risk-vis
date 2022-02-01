import React, { FC, useCallback } from 'react';
import { Box, Typography } from '@mui/material';

import { LAYERS } from '../config/layers';
import { CheckboxTree } from './checkbox-tree/CheckboxTree';
import { networkLayersConfig } from '../config/data/networks';
import _ from 'lodash';

interface NetworkControlProps {
  networkSelection: Record<string, boolean>;
  onNetworkSelection: (visUpdate: Record<string, boolean>) => void; //TODO change record key type to LayerName
}

const LayerLabel = ({ label, type, color }) => {
  return (
    <>
      <Typography>
        <Box component="span" className={type === 'line' ? 'dot line' : 'dot'} style={{ backgroundColor: color }}></Box>
        {label}
      </Typography>
    </>
  );
};
export const NetworkControl: FC<NetworkControlProps> = ({ networkSelection, onNetworkSelection }) => {
  const handleSelected = useCallback(
    (newSelected: string[]) => {
      onNetworkSelection({
        ..._.mapValues(networkSelection, () => false),
        ..._.fromPairs(newSelected.map((id) => [id, true])),
      });
    },
    [networkSelection, onNetworkSelection],
  );

  return (
    <Box mb={2}>
      <Typography variant="h6">Infrastructure Assets</Typography>
      <CheckboxTree
        nodes={networkLayersConfig}
        getLabel={(node) => (node.children ? node.label : <LayerLabel {...LAYERS[node.id]} label={node.label} />)}
        onSelected={handleSelected}
      />
    </Box>
  );
};
