import { FormControlLabel, Radio, RadioGroup } from '@mui/material';
import { useCallback } from 'react';
import { useRecoilState } from 'recoil';
import { boundaryLevelState } from './layers-state';

export const BoundaryLevelSelection = () => {
  const [boundaryLevel, setBoundaryLevel] = useRecoilState(boundaryLevelState);
  const handleBoundaryLevel = useCallback(
    (e, newBoundaryLevel) => {
      if (newBoundaryLevel != null) {
        setBoundaryLevel(newBoundaryLevel);
      }
    },
    [setBoundaryLevel],
  );
  return (
    <RadioGroup value={boundaryLevel} onChange={handleBoundaryLevel}>
      <FormControlLabel value="parish" label="Parishes" control={<Radio />} />
      {/* <FormControlLabel value="community" label="Communities" control={<Radio />} /> */}
      <FormControlLabel value="enumeration" label="Enumeration Districts" control={<Radio />} />
    </RadioGroup>
  );
};
