import { ComponentProps, FC, useCallback, useMemo } from 'react';
import { FormControl, FormLabel, makeStyles, Slider } from '@material-ui/core';

interface HazardsControlProps {
  fluvialReturnPeriod: number;
  setFluvialReturnPeriod: (years: number) => void;

  coastalReturnPeriod: number;
  setCoastalReturnPeriod: (years: number) => void;
}

interface GenericMark<T> {
  value: T;
  label: string;
}

type CustomSliderProps<T> = {
  marks: GenericMark<T>[];
  value: T;
  onChange: (newVal: T) => void;
} & Omit<ComponentProps<typeof Slider>, 'marks' | 'value' | 'onChange' | 'min' | 'max' | 'step' | 'scale'>;

const CustomNumberSlider: FC<CustomSliderProps<number>> = ({ marks, value, onChange, ...otherProps }) => {
  const integerMarks = useMemo(() => marks.map((m, idx) => ({ value: idx, label: m.label })), [marks]);
  const integerValue = useMemo(() => marks.findIndex((m) => m.value === value), [marks, value]);

  const handleIntegerChange = useCallback(
    (e, value: number) => {
      onChange(marks[value].value);
    },
    [marks, onChange],
  );

  return (
    <Slider
      marks={integerMarks}
      value={integerValue}
      onChange={handleIntegerChange}
      min={0}
      max={marks.length - 1}
      step={1}
      {...otherProps}
    />
  );
};

const fluvialMarks = [
  { value: 20, label: '20' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
  { value: 200, label: '200' },
  { value: 500, label: '500' },
  { value: 1500, label: '1500' },
];

const coastalMarks = [
  { value: 1, label: '1' },
  { value: 2, label: '2' },
  { value: 5, label: '5' },
  { value: 10, label: '10' },
  { value: 50, label: '50' },
  { value: 100, label: '100' },
];

const useStyles = makeStyles({
  sliderForm: {
    width: '100%',
  },
});

export const HazardsControl: FC<HazardsControlProps> = ({
  fluvialReturnPeriod,
  setFluvialReturnPeriod,
  coastalReturnPeriod,
  setCoastalReturnPeriod,
}) => {
  const classes = useStyles();

  return (
    <>
      <FormControl component="fieldset" className={classes.sliderForm}>
        <FormLabel component="legend">River Flooding Return Period</FormLabel>
        <CustomNumberSlider marks={fluvialMarks} value={fluvialReturnPeriod} onChange={setFluvialReturnPeriod} />
      </FormControl>
      <FormControl component="fieldset" className={classes.sliderForm}>
        <FormLabel component="legend">Coastal Flooding Return Period</FormLabel>
        <CustomNumberSlider marks={coastalMarks} value={coastalReturnPeriod} onChange={setCoastalReturnPeriod} />
      </FormControl>
    </>
  );
};
