import { Box, Button, Checkbox, FormControlLabel, Paper, ToggleButton, ToggleButtonGroup } from '@mui/material';
import { Layers as LayersIcon } from '@mui/icons-material';
import { useCallback, useState } from 'react';
import { useRecoilState } from 'recoil';

import { backgroundState, showLabelsState } from './layers-state';

export const MapLayerSelection = () => {
  const [showPopover, setShowPopover] = useState(false);

  const [background, setBackground] = useRecoilState(backgroundState);
  const [showLabels, setShowLabels] = useRecoilState(showLabelsState);

  const other = background === 'satellite' ? 'light' : 'satellite';

  const handleBackground = useCallback(
    (e, newBackground) => {
      // when user presses the currently active option
      if (newBackground == null) {
        // toggle to the other option
        setBackground(other);
      } else {
        setBackground(newBackground);
      }
    },
    [setBackground, other],
  );

  const buttonStyle = showPopover ? { borderTopRightRadius: 0, borderBottomRightRadius: 0 } : {};
  return (
    <Box
      onMouseEnter={() => setShowPopover(true)}
      onMouseLeave={() => setShowPopover(false)}
      style={{ display: 'flex' }}
    >
      <Button
        aria-label="Toggle map background"
        variant="contained"
        style={{
          paddingInline: 0,
          backgroundColor: 'white',
          color: 'black',
          minWidth: 'auto',
          width: '40px',
          height: '36px',
          ...buttonStyle,
        }}
      >
        <LayersIcon />
      </Button>
      {showPopover && (
        <Paper style={{ overflow: 'hidden', borderTopLeftRadius: 0 }}>
          <Box width={220}>
            <Box px={2} py="6px" height={37} display="flex" borderBottom="1px solid #ddd">
              <ToggleButtonGroup exclusive color="primary" value={background} onChange={handleBackground} fullWidth>
                <ToggleButton value="light">Map</ToggleButton>
                <ToggleButton value="satellite">Satellite</ToggleButton>
              </ToggleButtonGroup>
            </Box>
            <Box px={2}>
              <Box>
                <FormControlLabel
                  label="Show labels"
                  control={<Checkbox checked={showLabels} onChange={(e, checked) => setShowLabels(checked)} />}
                />
              </Box>
            </Box>
          </Box>
        </Paper>
      )}
    </Box>
  );
};
