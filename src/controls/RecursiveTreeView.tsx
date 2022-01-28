import React, { useCallback, useMemo } from 'react';
import { TreeView } from '@material-ui/lab';
import { ExpandMore as ExpandMoreIcon, ChevronRight as ChevronRightIcon } from '@material-ui/icons';
import TreeItem from '@material-ui/lab/TreeItem';
import { Checkbox, FormControlLabel } from '@material-ui/core';

export type TreeNode<T> = {
  id: string;
  children?: TreeNode<T>[];
} & T;

function getChildrenById<T>(node: TreeNode<T>, id: string) {
  let array: string[] = [];

  function getAllChildren(nodes: TreeNode<T> | null) {
    if (nodes === null) return [];
    array.push(nodes.id);
    if (Array.isArray(nodes.children)) {
      nodes.children.forEach((node) => {
        array = [...array, ...getAllChildren(node)];
        array = array.filter((v, i) => array.indexOf(v) === i);
      });
    }
    return array;
  }

  function getNodeById<T>(nodes: TreeNode<T>, id: string) {
    if (nodes.id === id) {
      return nodes;
    } else if (Array.isArray(nodes.children)) {
      let result = null;
      nodes.children.forEach((node) => {
        if (getNodeById(node, id)) {
          result = getNodeById(node, id);
        }
      });
      return result;
    }

    return null;
  }

  return getAllChildren(getNodeById(node, id));
}

function RecursiveTreeItem<T>({
  root,
  selected,
  handleChange,
  getLabel,
}: {
  root: TreeNode<T>;
  selected: string[];
  handleChange: (checked: boolean, node: TreeNode<T>) => void;
  getLabel: (node: TreeNode<T>) => any;
}) {
  return (
    <TreeItem
      key={root.id}
      nodeId={root.id}
      label={
        <FormControlLabel
          key={root.id}
          label={getLabel(root)}
          style={{ pointerEvents: 'none' }}
          control={
            <Checkbox
              checked={selected.some((item) => item === root.id)}
              onChange={(event) => handleChange(event.currentTarget.checked, root)}
              onClick={(e) => e.stopPropagation()}
              style={{ pointerEvents: 'auto' }}
            />
          }
        ></FormControlLabel>
      }
    >
      {Array.isArray(root.children)
        ? root.children.map((node) => (
            <RecursiveTreeItem
              root={node}
              selected={selected}
              handleChange={handleChange}
              getLabel={getLabel}
            ></RecursiveTreeItem>
          ))
        : null}
    </TreeItem>
  );
}

function buildConfig<T>(nodes: TreeNode<T>[]) {
  function dfs(node: TreeNode<T>, action: (node: TreeNode<T>) => void) {
    action(node);
    if (Array.isArray(node.children)) {
      node.children.forEach((child) => dfs(child, action));
    }
  }

  const config = {};

  nodes.forEach((node) => {
    dfs(node, (node) => {
      config[node.id] = {
        ...node,
      };
    });
  });
  return config;
}

export default function RecursiveTreeView<T>({
  nodes,
  getLabel,
  onSelected,
}: {
  nodes: TreeNode<T>[];
  getLabel: (node: TreeNode<T>) => any;
  onSelected: (selected: string[]) => void;
}) {
  const config = useMemo(() => buildConfig(nodes), [nodes]);
  const [selected, setSelected] = React.useState<string[]>([]);

  const handleSelected = useCallback(
    (newSelected: string[]) => {
      setSelected(newSelected);

      const selectedLeafNodes = newSelected.filter((x) => !config[x].children);

      onSelected(selectedLeafNodes);
    },
    [config, onSelected],
  );

  function handleChange(checked: boolean, nodes: TreeNode<T>) {
    const allNode: string[] = getChildrenById(nodes, nodes.id);
    let array = checked ? [...selected, ...allNode] : selected.filter((value) => !allNode.includes(value));

    array = array.filter((v, i) => array.indexOf(v) === i);

    handleSelected(array);
  }

  return (
    <>
      <TreeView
        defaultCollapseIcon={<ExpandMoreIcon />}
        defaultExpanded={['0', '3', '4']}
        defaultExpandIcon={<ChevronRightIcon />}
        style={{
          flexGrow: 1,
          maxWidth: 400,
        }}
      >
        {nodes.map((node) => (
          <RecursiveTreeItem
            root={node}
            selected={selected}
            handleChange={handleChange}
            getLabel={getLabel}
          ></RecursiveTreeItem>
        ))}
      </TreeView>
    </>
  );
}
