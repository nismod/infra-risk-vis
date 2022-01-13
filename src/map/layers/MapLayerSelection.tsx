import { Button, ButtonGroup, Typography } from '@material-ui/core';
import { Layers as LayersIcon } from '@material-ui/icons';
import { useState } from 'react';
import { useRecoilState } from 'recoil';

import { backgroundState } from './layers-state';

const config = {
  satellite: {
    label: 'Satellite',
  },
  light: {
    label: 'Map',
  },
};

export const MapLayerSelection = () => {
  const [background, setBackground] = useRecoilState(backgroundState);
  const [showHint, setShowHint] = useState(false);

  const other = background === 'satellite' ? 'light' : 'satellite';

  // const currentConfig = config[background];
  return (
    <ButtonGroup>
      <Button
        aria-label="Toggle map background"
        onClick={() => setBackground(other)}
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
  );
};
