import { useRecoilValue } from 'recoil';

import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { Paper, Box } from '@mui/material';
import { RegionDetailsContent } from './RegionDetailsContent';

export const RegionDetails = () => {
  const selectedRegion = useRecoilValue(selectionState('regions'));

  if (!selectedRegion) return null;

  return (
    <Paper>
      <Box p={3} maxHeight="calc(100vh - 125px)" style={{ overflowY: 'scroll' }}>
        <RegionDetailsContent selectedRegion={selectedRegion} />
      </Box>
    </Paper>
  );
};
