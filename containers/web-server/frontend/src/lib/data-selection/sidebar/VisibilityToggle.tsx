import { Visibility, VisibilityOff } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { FC } from 'react';

export interface VisibilityToggleProps {
  visibility: boolean;
  onVisibility: (x: boolean) => void;
  labelShow?: string;
  labelHide?: string;
}

export const VisibilityToggle: FC<VisibilityToggleProps> = ({
  visibility,
  onVisibility,
  labelShow = 'Show layer',
  labelHide = 'Hide layer',
}) => {
  return (
    <IconButton
      size="small"
      title={visibility ? labelHide : labelShow}
      onClick={(e) => {
        onVisibility(!visibility);
        e.stopPropagation();
      }}
    >
      {visibility ? <Visibility color="action" /> : <VisibilityOff color="disabled" />}
    </IconButton>
  );
};
