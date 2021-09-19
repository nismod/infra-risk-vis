import { ChangeEvent, ComponentProps, FC, useCallback, useMemo } from 'react';
import { Slider } from '@material-ui/core';

export interface GenericMark<T> {
  value: T;
  label: string;
}

type CustomSliderProps<T> = {
  marks: GenericMark<T>[];
  value: T;
  onChange: (newVal: T) => void;
} & Omit<ComponentProps<typeof Slider>, 'marks' | 'value' | 'onChange' | 'min' | 'max' | 'step' | 'scale'>;

export const CustomNumberSlider: FC<CustomSliderProps<number>> = ({ marks, value, onChange, ...otherProps }) => {
  const integerMarks = useMemo(() => marks.map((m, idx) => ({ value: idx, label: m.label })), [marks]);
  const integerValue = useMemo(() => marks.findIndex((m) => m.value === value), [marks, value]);

  const handleIntegerChange = useCallback(
    (e: ChangeEvent<{}>, value: number) => {
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
