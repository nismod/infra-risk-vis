import { Box, Paper } from '@material-ui/core';

export const MapLegend = ({ children }) => {
  return (
    <Paper>
      <Box p={1}>{children}</Box>
    </Paper>
  );
};
