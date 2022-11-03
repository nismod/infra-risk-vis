import _ from 'lodash';
import { atom, selector } from 'recoil';

import { CheckboxTreeState, buildTreeConfig } from '@/lib/controls/checkbox-tree/CheckboxTree';

import { NETWORK_LAYERS_HIERARCHY } from '@/config/networks/hierarchy';
import { NetworkLayerType } from '@/config/networks/metadata';

export const networkTreeExpandedState = atom<string[]>({
  key: 'networkTreeExpandedState',
  default: [],
});

export const networkTreeConfig = buildTreeConfig(NETWORK_LAYERS_HIERARCHY);

export const networkTreeCheckboxState = atom<CheckboxTreeState>({
  key: 'networkTreeSelectionState',
  default: {
    checked: _.mapValues(networkTreeConfig.nodes, () => false),
    indeterminate: _.mapValues(networkTreeConfig.nodes, () => false),
  },
});

export const networkSelectionState = selector<NetworkLayerType[]>({
  key: 'networkSelectionState',
  get: ({ get }) => {
    const checkboxState = get(networkTreeCheckboxState);

    return Object.keys(checkboxState.checked).filter(
      (id) => checkboxState.checked[id] && !networkTreeConfig.nodes[id].children,
    ) as NetworkLayerType[];
  },
});
