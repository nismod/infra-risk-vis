import { BoxProps, Paper } from '@mui/material';
import { Box } from '@mui/system';
import { FC } from 'react';

export const SidePanel: FC<BoxProps> = ({ children, ...props }) => {
  return (
    <Paper>
      <Box px={3} py={2} {...props}>
        {children}
      </Box>
    </Paper>
  );
};
