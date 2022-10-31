import { TreeItem } from '@mui/lab';
import { Checkbox, FormControlLabel } from '@mui/material';

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
  getLabel: (node: TreeNode<T>) => any;
  disableCheck?: boolean;
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
              checked={checkboxState.indeterminate[root.id] || checkboxState.checked[root.id]}
              indeterminate={checkboxState.indeterminate[root.id]}
              onChange={(event) => handleChange(event.currentTarget.checked, root)}
              onClick={(e) => e.stopPropagation()}
              style={{ pointerEvents: 'auto' }}
              disabled={disableCheck}
            />
          }
        ></FormControlLabel>
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
