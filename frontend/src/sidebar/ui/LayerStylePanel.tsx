import { Paper } from '@mui/material';
import { Box } from '@mui/system';
import { FC } from 'react';

export const LayerStylePanel: FC<{}> = ({ children }) => {
  return (
    <Box mt={2}>
      <Paper>
        <Box p={2}>{children}</Box>
      </Paper>
    </Box>
  );
};
