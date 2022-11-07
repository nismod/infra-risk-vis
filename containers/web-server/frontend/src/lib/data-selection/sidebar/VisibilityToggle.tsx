import { Visibility, VisibilityOff } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { FC } from 'react';

export interface VisibilityToggleProps {
  visibility: boolean;
  onVisibility: (x: boolean) => void;
}

export const VisibilityToggle: FC<VisibilityToggleProps> = ({ visibility, onVisibility }) => {
  return (
    <IconButton
      size="small"
      title={visibility ? 'Hide layer' : 'Show layer'}
      onClick={(e) => {
        onVisibility(!visibility);
        e.stopPropagation();
      }}
    >
      {visibility ? <Visibility color="action" /> : <VisibilityOff color="disabled" />}
    </IconButton>
  );
};
