import { TreeItem } from '@mui/lab';
import { Box, Checkbox } from '@mui/material';

import { CheckboxTreeState } from './CheckboxTree';
import { TreeNode } from './tree-node';

export function CheckboxTreeItem<T>({
  root,
  handleChange,
  checkboxState,
  getLabel,
  disableCheck = false,
}: {
  root: TreeNode<T>;
  handleChange: (checked: boolean, node: TreeNode<T>) => void;
  checkboxState: CheckboxTreeState;
  getLabel: (node: TreeNode<T>, checked: boolean) => any;
  disableCheck?: boolean;
}) {
  const indeterminate = checkboxState.indeterminate[root.id];
  const checked = indeterminate || checkboxState.checked[root.id];
  return (
    <TreeItem
      key={root.id}
      nodeId={root.id}
      label={
        <Box display="flex" alignItems="center" width="100%">
          <Checkbox
            checked={checked}
            indeterminate={indeterminate}
            onChange={(event) => handleChange(event.currentTarget.checked, root)}
            onClick={(e) => e.stopPropagation()}
            disabled={disableCheck}
          />
          <Box flexGrow={1}>{getLabel(root, checked)}</Box>
        </Box>
      }
    >
      {root.children?.map((node) => (
        <CheckboxTreeItem
          key={node.id}
          root={node}
          handleChange={handleChange}
          checkboxState={checkboxState}
          getLabel={getLabel}
          disableCheck={disableCheck}
        ></CheckboxTreeItem>
      ))}
    </TreeItem>
  );
}
