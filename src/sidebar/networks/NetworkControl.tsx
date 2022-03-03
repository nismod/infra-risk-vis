import _ from 'lodash';
import { FC, useCallback } from 'react';
import { Box, Typography } from '@mui/material';

import { CheckboxTree } from 'lib/controls/checkbox-tree/CheckboxTree';

import { NETWORK_LAYERS_HIERARCHY } from 'config/networks/hierarchy';
import { useRecoilState } from 'recoil';
import { networkSelectionState } from 'state/network-selection';
import { NETWORKS_METADATA } from 'config/networks/metadata';

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
    <CheckboxTree
      nodes={NETWORK_LAYERS_HIERARCHY}
      getLabel={(node) =>
        node.children ? node.label : <LayerLabel {...NETWORKS_METADATA[node.id]} label={node.label} />
      }
      onSelected={handleSelected}
    />
  );
};
