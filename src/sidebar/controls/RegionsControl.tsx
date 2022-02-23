import { Checkbox, FormControlLabel } from '@mui/material';
import { BoundaryLevelSelection } from 'map/layers/BoundaryLevelSelection';
import { useRecoilState, useRecoilValue } from 'recoil';
import { showPopulationState } from 'state/population';

const PopulationToggle = () => {
  const [showPopulation, setShowPopulation] = useRecoilState(showPopulationState);

  return (
    <FormControlLabel
      label="Show Population Map"
      control={<Checkbox checked={showPopulation} onChange={(e, checked) => setShowPopulation(checked)} />}
    />
  );
};

export const RegionsControl = () => {
  const showPopulation = useRecoilValue(showPopulationState);

  return (
    <>
      <PopulationToggle />
      {showPopulation && <BoundaryLevelSelection />}
    </>
  );
};
