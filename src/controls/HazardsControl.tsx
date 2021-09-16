import { FC, useCallback, useMemo } from 'react';
import { FormControl, FormLabel, makeStyles, Slider } from '@material-ui/core';

interface HazardsControlProps {
  returnPeriodYears: number;
  setReturnPeriodYears: (years: number) => void;
}

interface GenericMark<T> {
  value: T;
  label: string;
}

function useCustomSlider<T>(marks: GenericMark<T>[], value: T, onChange: (newVal: T) => void) {
  const integerMarks = useMemo(() => marks.map((m, idx) => ({ value: idx, label: m.label })), [marks]);
  const integerValue = useMemo(() => marks.findIndex((m) => m.value === value), [marks, value]);

  const handleIntegerChange = useCallback(
    (e, value: number) => {
      onChange(marks[value].value);
    },
    [marks, onChange],
  );

  return {
    integerMarks,
    integerValue,
    handleIntegerChange,
  };
}

const marks = [
  { value: 20, label: '20' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
  { value: 200, label: '200' },
  { value: 500, label: '500' },
  { value: 1500, label: '1500' },
];

const useStyles = makeStyles({
  sliderForm: {
    width: '100%',
  },
});

export const HazardsControl: FC<HazardsControlProps> = ({ returnPeriodYears, setReturnPeriodYears }) => {
  const { integerMarks, integerValue, handleIntegerChange } = useCustomSlider(
    marks,
    returnPeriodYears,
    setReturnPeriodYears,
  );

  const classes = useStyles();

  return (
    <>
      <FormControl component="fieldset" className={classes.sliderForm}>
        <FormLabel component="legend">River Flooding Return Period</FormLabel>
        <Slider
          marks={integerMarks}
          value={integerValue}
          onChange={handleIntegerChange}
          min={0}
          max={marks.length - 1}
          step={1}
        />
      </FormControl>
    </>
  );
};
