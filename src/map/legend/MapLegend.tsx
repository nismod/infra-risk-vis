import { Box, Paper } from '@mui/material';

export const MapLegend = ({ children }) => {
  return (
    <Paper>
      <Box p={1}>{children}</Box>
    </Paper>
  );
};
