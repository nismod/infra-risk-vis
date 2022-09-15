import { ComponentProps, FC, useCallback, useMemo } from 'react';
import { Slider } from '@mui/material';

export interface GenericMark<T> {
  value: T;
  label: string;
}

type CustomSliderProps<T> = {
  marks: T[]; //GenericMark<T>[];
  value: T;
  onChange: (newVal: T) => void;
  showMarkLabelsFor?: T[];
} & Omit<ComponentProps<typeof Slider>, 'marks' | 'value' | 'onChange' | 'min' | 'max' | 'step' | 'scale'>;

export const CustomNumberSlider: FC<CustomSliderProps<number>> = ({
  marks,
  value,
  onChange,
  showMarkLabelsFor,
  ...otherProps
}) => {
  const showLabelLookup = useMemo(
    () => showMarkLabelsFor && Object.fromEntries(showMarkLabelsFor.map((ml) => [ml, true])),
    [showMarkLabelsFor],
  );
  const integerMarks = useMemo(() => {
    return marks.map((m, idx) => ({
      value: idx,
      label: !showLabelLookup || showLabelLookup[m] ? m.toString() : null,
    }));
  }, [marks, showLabelLookup]);
  const integerValue = useMemo(() => marks.findIndex((m) => m === value), [marks, value]);

  const handleIntegerChange = useCallback(
    (e: Event, value: number) => {
      onChange(marks[value]);
    },
    [marks, onChange],
  );

  const valueLabelFunction = useCallback((value) => marks[value].toString(), [marks]);

  return (
    <Slider
      marks={integerMarks}
      value={integerValue}
      onChange={handleIntegerChange}
      min={0}
      max={marks.length - 1}
      step={1}
      valueLabelFormat={valueLabelFunction}
      getAriaValueText={valueLabelFunction}
      {...otherProps}
    />
  );
};
