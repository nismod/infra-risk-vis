import { Alert } from '@mui/material';
import { Box } from '@mui/system';
import { FC } from 'react';
import { useRecoilState, useRecoilValue } from 'recoil';

import { CheckboxTree } from '@/lib/controls/checkbox-tree/CheckboxTree';

import { NETWORK_LAYERS_HIERARCHY } from '@/config/networks/hierarchy';
import { NETWORKS_METADATA } from '@/config/networks/metadata';
import { LayerLabel } from '@/sidebar/ui/LayerLabel';
import { showAdaptationsState } from '@/state/data-selection/networks/adaptations';
import {
  networkTreeCheckboxState,
  networkTreeConfig,
  networkTreeExpandedState,
} from '@/state/data-selection/networks/network-selection';

export const NetworkControl: FC<{}> = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(networkTreeCheckboxState);
  const [expanded, setExpanded] = useRecoilState(networkTreeExpandedState);

  const showAdaptations = useRecoilValue(showAdaptationsState);
  const disableCheck = showAdaptations;

  return (
    <>
      {showAdaptations ? (
        <Box my={1}>
          <Alert severity="info">
            Infrastructure layers are currently following the Adaptation Options selection
          </Alert>
        </Box>
      ) : null}
      <CheckboxTree
        nodes={NETWORK_LAYERS_HIERARCHY}
        config={networkTreeConfig}
        getLabel={(node) =>
          node.children ? (
            node.label
          ) : (
            <LayerLabel {...NETWORKS_METADATA[node.id]} label={node.label} />
          )
        }
        checkboxState={checkboxState}
        onCheckboxState={setCheckboxState}
        expanded={expanded}
        onExpanded={setExpanded}
        disableCheck={disableCheck}
      />
    </>
  );
};
