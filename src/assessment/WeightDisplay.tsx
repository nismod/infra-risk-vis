import { Slider } from './ValueDisplay';
import { numRound } from 'lib/helpers';

export const WeightDisplay = ({ value }: { value: number }) => {
  return (
    <Slider
      disabled
      min={0}
      max={1}
      marks={[
        { value: 0, label: '0' },
        { value: 1, label: '+' },
      ]}
      track={false}
      value={numRound(value)}
      valueLabelDisplay="on"
    />
  );
};
