import { Box, Button, ButtonGroup, Typography } from '@material-ui/core';
import { Layers as LayersIcon } from '@material-ui/icons';
import { useState } from 'react';

const config = {
  satellite: {
    label: 'Satellite',
  },
  light: {
    label: 'Map',
  },
};

export const MapLayerSelection = ({ background, onBackground }) => {
  const [showHint, setShowHint] = useState(false);

  const other = background === 'satellite' ? 'light' : 'satellite';

  // const currentConfig = config[background];
  return (
    <Box position="absolute" top={0} left={0} m={1} zIndex={1000}>
      <ButtonGroup>
        <Button
          aria-label="Toggle map background"
          onClick={() => onBackground(other)}
          onMouseEnter={() => setShowHint(true)}
          onMouseLeave={() => setShowHint(false)}
          variant="contained"
          style={{ paddingInline: 0, backgroundColor: 'white' }}
        >
          <LayersIcon />
        </Button>
        {showHint && (
          <Button variant="contained" color="default" style={{ textTransform: 'none' }}>
            <Typography>Switch to {config[other].label} background</Typography>
          </Button>
        )}
      </ButtonGroup>
    </Box>
  );
};
