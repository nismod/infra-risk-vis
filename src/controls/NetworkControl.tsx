import _ from 'lodash';
import { FC, useCallback } from 'react';
import { Box, Typography } from '@mui/material';

import { CheckboxTree } from 'lib/controls/checkbox-tree/CheckboxTree';

import { LAYERS } from '../config/layers';
import { networkLayersConfig } from '../config/data/networks';
import { useRecoilState } from 'recoil';
import { networkSelectionState } from 'state/network-selection';

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
export const NetworkControl: FC<{}> = () => {
  const [networkSelection, setNetworkSelection] = useRecoilState(networkSelectionState);

  const handleSelected = useCallback(
    (newSelected: string[]) => {
      setNetworkSelection({
        ..._.mapValues(networkSelection, () => false),
        ..._.fromPairs(newSelected.map((id) => [id, true])),
      });
    },
    [networkSelection, setNetworkSelection],
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
