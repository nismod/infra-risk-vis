import { Close } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { useResetRecoilState } from 'recoil';

import { selectionState } from '@/lib/data-map/interactions/interaction-state';

export const DeselectButton = ({ interactionGroup, title }) => {
  const clearSelection = useResetRecoilState(selectionState(interactionGroup));

  return (
    <IconButton onClick={() => clearSelection()} title={title}>
      <Close />
    </IconButton>
  );
};
