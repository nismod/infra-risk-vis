import { useRecoilValue, useResetRecoilState } from 'recoil';

import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { Paper, Box, IconButton } from '@mui/material';
import { RegionDetailsContent } from './RegionDetailsContent';
import { Close } from '@mui/icons-material';

export const RegionDetails = () => {
  const selectedRegion = useRecoilValue(selectionState('regions'));
  const clearSelectedRegion = useResetRecoilState(selectionState('regions'));

  if (!selectedRegion) return null;

  return (
    <Paper>
      <Box p={3} maxHeight="30vh" style={{ overflowY: 'scroll' }} position="relative">
        <Box position="absolute" top={0} right={0} p={2}>
          <IconButton onClick={() => clearSelectedRegion()} title="Deselect region">
            <Close />
          </IconButton>
        </Box>
        <RegionDetailsContent selectedRegion={selectedRegion} />
      </Box>
    </Paper>
  );
};
