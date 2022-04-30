import { useRecoilValue, useResetRecoilState } from 'recoil';

import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { Box, IconButton } from '@mui/material';
import { RegionDetailsContent } from './RegionDetailsContent';
import { Close } from '@mui/icons-material';
import { SidePanel } from 'details/SidePanel';

export const RegionDetails = () => {
  const selectedRegion = useRecoilValue(selectionState('regions'));
  const clearSelectedRegion = useResetRecoilState(selectionState('regions'));

  if (!selectedRegion) return null;

  return (
    <SidePanel position="relative">
      <Box position="absolute" top={0} right={0} p={2}>
        <IconButton onClick={() => clearSelectedRegion()} title="Deselect region">
          <Close />
        </IconButton>
      </Box>
      <RegionDetailsContent selectedRegion={selectedRegion} />
    </SidePanel>
  );
};
