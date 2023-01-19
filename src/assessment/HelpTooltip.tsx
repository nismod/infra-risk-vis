import { Help } from '@mui/icons-material';
import { IconButton, Tooltip } from '@mui/material';

export function HelpTooltip(text: string) {
  return (
    <Tooltip title={text}>
      <IconButton size="small">
        <Help fontSize="inherit" />
      </IconButton>
    </Tooltip>
  );
}
