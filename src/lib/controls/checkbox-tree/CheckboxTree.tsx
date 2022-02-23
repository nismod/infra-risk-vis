import React, { useMemo } from 'react';
import _ from 'lodash';
import { TreeView } from '@mui/lab';
import { ExpandMore as ExpandMoreIcon, ChevronRight as ChevronRightIcon } from '@mui/icons-material';
import { useImmer } from 'use-immer';

import { useChangeEffect } from 'lib/hooks/use-change-effect';

import { dfs, getDescendants, TreeNode } from './tree-node';
import { CheckboxTreeItem } from './CheckboxTreeItem';

export interface CheckboxTreeConfig<T> {
  roots: TreeNode<T>[];
  nodes: {
    [nodeId: string]: TreeNode<T> & {
      descendants: string[];
    };
  };
}

function buildConfig<T>(nodes: TreeNode<T>[]): CheckboxTreeConfig<T> {
  const config: CheckboxTreeConfig<T> = {
    roots: nodes,
    nodes: {},
  };

  nodes.forEach((node) => {
    dfs(node, (node) => {
      config.nodes[node.id] = {
        ...node,
        descendants: getDescendants(node),
      };
    });
  });
  return config;
}

export interface CheckboxTreeState {
  checked: { [nodeId: string]: boolean };
  indeterminate: { [nodeId: string]: boolean };
}

function initState<T>(config: CheckboxTreeConfig<T>): CheckboxTreeState {
  return {
    checked: _.mapValues(config.nodes, () => false),
    indeterminate: _.mapValues(config.nodes, () => false),
  };
}

function recalculateCheckboxStates<T>(state: CheckboxTreeState, config: CheckboxTreeConfig<T>): CheckboxTreeState {
  for (const root of config.roots) {
    // traverse each root tree in post-order to recalculate state starting from leaf nodes
    dfs(
      root,
      (node) => {
        const nodeChildren = config.nodes[node.id].children;
        if (nodeChildren) {
          const checked = nodeChildren.every((child) => state.checked[child.id]);
          const indeterminate =
            !checked && nodeChildren.some((child) => state.checked[child.id] || state.indeterminate[child.id]);
          state.checked[node.id] = checked;
          state.indeterminate[node.id] = indeterminate;
        }
      },
      false,
      'post',
    );
  }

  return state;
}

export function CheckboxTree<T>({
  nodes,
  getLabel,
  onSelected,
}: {
  nodes: TreeNode<T>[];
  getLabel: (node: TreeNode<T>) => any;
  onSelected: (selected: string[]) => void;
}) {
  const config = useMemo(() => buildConfig(nodes), [nodes]);
  const [checkboxState, setCheckboxState] = useImmer(() => initState(config));

  useChangeEffect(
    () => {
      const selected = Object.keys(checkboxState.checked).filter(
        (id) => checkboxState.checked[id] && !config.nodes[id].children,
      );
      onSelected(selected);
    },
    [checkboxState, config, onSelected],
    [checkboxState],
  );

  function handleChange(checked: boolean, node: TreeNode<T>) {
    const descendants: string[] = config.nodes[node.id].descendants;
    setCheckboxState((draft) => {
      draft.checked[node.id] = checked;
      descendants.forEach((n) => (draft.checked[n] = checked));

      return recalculateCheckboxStates(draft, config);
    });
  }

  return (
    <>
      <TreeView defaultCollapseIcon={<ExpandMoreIcon />} defaultExpandIcon={<ChevronRightIcon />}>
        {nodes.map((node) => (
          <CheckboxTreeItem
            key={node.id}
            root={node}
            checkboxState={checkboxState}
            handleChange={handleChange}
            getLabel={getLabel}
          />
        ))}
      </TreeView>
    </>
  );
}
