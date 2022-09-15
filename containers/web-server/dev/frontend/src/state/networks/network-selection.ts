import { NETWORK_LAYERS_HIERARCHY } from 'config/networks/hierarchy';
import { buildTreeConfig, CheckboxTreeState } from 'lib/controls/checkbox-tree/CheckboxTree';
import _ from 'lodash';
import { atom, selector } from 'recoil';

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

export const networkSelectionState = selector<string[]>({
  key: 'networkSelectionState',
  get: ({ get }) => {
    const checkboxState = get(networkTreeCheckboxState);

    return Object.keys(checkboxState.checked).filter(
      (id) => checkboxState.checked[id] && !networkTreeConfig.nodes[id].children,
    );
  },
});
