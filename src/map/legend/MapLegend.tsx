import { Box, Paper } from '@material-ui/core';

export const MapLegend = ({ children }) => {
  return (
    <Box position="absolute" bottom={0} left={0} m={1} zIndex={1000}>
      <Paper>
        <Box p={1}>{children}</Box>
      </Paper>
    </Box>
  );
};
