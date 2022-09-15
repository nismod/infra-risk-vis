import { Visibility, VisibilityOff } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { useRecoilState } from 'recoil';
import { sectionVisibilityState } from 'state/sections';

export const VisibilityToggle = ({ id }) => {
  const [visibility, setVisibility] = useRecoilState(sectionVisibilityState(id));

  return (
    <IconButton
      title={visibility ? 'Hide layer' : 'Show layer'}
      onClick={(e) => {
        setVisibility((visibility) => !visibility);
        e.stopPropagation();
      }}
    >
      {visibility ? <Visibility /> : <VisibilityOff />}
    </IconButton>
  );
};
