import { Checkbox, FormControlLabel } from '@material-ui/core';
import { TreeItem } from '@material-ui/lab';
import { CheckboxTreeState } from './CheckboxTree';
import { TreeNode } from './tree-node';

export function CheckboxTreeItem<T>({
  root,
  handleChange,
  checkboxState,
  getLabel,
}: {
  root: TreeNode<T>;
  handleChange: (checked: boolean, node: TreeNode<T>) => void;
  checkboxState: CheckboxTreeState;
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
              checked={checkboxState.indeterminate[root.id] || checkboxState.checked[root.id]}
              indeterminate={checkboxState.indeterminate[root.id]}
              onChange={(event) => handleChange(event.currentTarget.checked, root)}
              onClick={(e) => e.stopPropagation()}
              style={{ pointerEvents: 'auto' }}
            />
          }
        ></FormControlLabel>
      }
    >
      {root.children?.map((node) => (
        <CheckboxTreeItem
          root={node}
          handleChange={handleChange}
          checkboxState={checkboxState}
          getLabel={getLabel}
        ></CheckboxTreeItem>
      ))}
    </TreeItem>
  );
}
