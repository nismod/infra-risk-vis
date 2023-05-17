import { useContext } from 'react';

import { ViewStateContext } from 'lib/data-map/DeckMap';
import { Box, Paper } from '@mui/material';

export const ViewStateDebug = () => {
  const { viewState } = useContext(ViewStateContext);

  return (
    <Paper sx={{ my: 1 }}>
      <Box maxWidth="200px" px={2}>
        <p>Zoom: {viewState.zoom.toFixed(2)}</p>
        <p>Lat: {viewState.latitude.toFixed(4)}</p>
        <p>Lon: {viewState.longitude.toFixed(4)}</p>
      </Box>
    </Paper>
  );
};
