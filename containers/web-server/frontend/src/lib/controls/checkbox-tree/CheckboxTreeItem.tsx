import { TreeItem } from '@mui/lab';
import { Checkbox, FormControlLabel } from '@mui/material';
import { MouseEventHandler, useCallback } from 'react';

import { CheckboxTreeState } from './CheckboxTree';
import { TreeNode } from './tree-node';

export function CheckboxTreeItem<T>({
  root,
  handleChange,
  checkboxState,
  getLabel,
  disableCheck = false,
  toggleOnLeafClick = true,
}: {
  root: TreeNode<T>;
  handleChange: (checked: boolean, node: TreeNode<T>) => void;
  checkboxState: CheckboxTreeState;
  getLabel: (node: TreeNode<T>) => any;
  disableCheck?: boolean;
  /**
   * When clicking on a leaf item, should the checkbox be toggled?
   * @defaultValue true
   */
  toggleOnLeafClick?: boolean;
}) {
  // set checked to true if indeterminate, so that clicking on an indeterminate node toggles it off, not on
  const effectiveChecked = checkboxState.indeterminate[root.id] || checkboxState.checked[root.id];

  const isLeaf = !root.children || root.children.length === 0;
  const handleLeafClick = useCallback<MouseEventHandler<HTMLLIElement>>(
    (e) => {
      handleChange(!effectiveChecked, root);
    },
    [effectiveChecked, handleChange, root],
  );

  return (
    <TreeItem
      key={root.id}
      nodeId={root.id}
      onClick={toggleOnLeafClick && isLeaf ? handleLeafClick : undefined}
      label={
        <FormControlLabel
          key={root.id}
          label={getLabel(root)}
          style={{ pointerEvents: 'none' }}
          control={
            <Checkbox
              checked={effectiveChecked}
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
          toggleOnLeafClick={toggleOnLeafClick}
        ></CheckboxTreeItem>
      ))}
    </TreeItem>
  );
}
