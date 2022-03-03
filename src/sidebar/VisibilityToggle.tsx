import { Visibility, VisibilityOff } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { useRecoilState } from 'recoil';
import { showLayerState } from 'state/show-layers';

export const VisibilityToggle = ({ id }) => {
  const [visibility, setVisibility] = useRecoilState(showLayerState(id));

  return (
    <IconButton
      onClick={(e) => {
        setVisibility((visibility) => !visibility);
        e.stopPropagation();
      }}
    >
      {visibility ? <Visibility /> : <VisibilityOff />}
    </IconButton>
  );
};
