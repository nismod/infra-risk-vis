import { Box, Button, Checkbox, FormControlLabel, Paper, Typography } from '@material-ui/core';
import { Layers as LayersIcon } from '@material-ui/icons';
import { useState } from 'react';
import { useRecoilState } from 'recoil';

import { backgroundState, showLabelsState, showRegionsState } from './layers-state';

const config = {
  satellite: {
    label: 'Satellite',
  },
  light: {
    label: 'Map',
  },
};

export const MapLayerSelection = () => {
  const [showPopover, setShowPopover] = useState(false);

  const [background, setBackground] = useRecoilState(backgroundState);
  const [showLabels, setShowLabels] = useRecoilState(showLabelsState);
  const [showRegions, setShowRegions] = useRecoilState(showRegionsState);

  const other = background === 'satellite' ? 'light' : 'satellite';

  const buttonStyle = showPopover ? { borderTopRightRadius: 0, borderBottomRightRadius: 0 } : {};
  return (
    <Box
      onMouseEnter={() => setShowPopover(true)}
      onMouseLeave={() => setShowPopover(false)}
      style={{ display: 'flex', pointerEvents: 'none' }}
    >
      <Button
        aria-label="Toggle map background"
        onClick={() => setBackground(other)}
        variant="contained"
        style={{
          paddingInline: 0,
          backgroundColor: 'white',
          minWidth: 'auto',
          width: '40px',
          height: '36px',
          pointerEvents: 'auto',
          ...buttonStyle,
        }}
      >
        <LayersIcon />
      </Button>
      {showPopover && (
        <Paper style={{ overflow: 'hidden', pointerEvents: 'auto', borderTopLeftRadius: 0 }}>
          <Box px={2} py="6px" height={36} bgcolor="#ddd" display="flex">
            <Typography>Switch to {config[other].label} background</Typography>
          </Box>
          <Box px={2}>
            <Box>
              <FormControlLabel
                label="Show labels"
                control={<Checkbox checked={showLabels} onChange={(e, checked) => setShowLabels(checked)} />}
              />
            </Box>
            <Box>
              <FormControlLabel
                label="Show boundaries"
                control={<Checkbox checked={showRegions} onChange={(e, checked) => setShowRegions(checked)} />}
              />
            </Box>
          </Box>
        </Paper>
      )}
    </Box>
  );
};
