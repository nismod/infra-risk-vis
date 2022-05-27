import { Box } from '@mui/material';
import { FC } from 'react';

export const SidebarPanelSection: FC<{
  variant?: 'standard' | 'style';
}> = ({ children, variant = 'standard' }) => {
  const bgcolor = variant === 'style' ? '#eee' : undefined;
  return (
    <Box p={2} bgcolor={bgcolor}>
      {children}
    </Box>
  );
};
