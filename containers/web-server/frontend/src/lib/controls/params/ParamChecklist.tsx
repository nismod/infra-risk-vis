import { Button, Checkbox, FormControlLabel, FormGroup, FormLabel, Stack } from '@mui/material';
import { Box } from '@mui/system';
import { PropsWithChildren, ReactElement, useMemo } from 'react';

import { fromKeys } from '@/lib/helpers';

import { ValueLabel, getValueLabel } from './value-label';

interface ParamChecklistProps<K extends string> {
  title: string;
  options: K[] | ValueLabel<K>[];
  checklistState: { [key in K]: boolean };
  onChecklistState: (state: { [key in K]: boolean }) => void;
  renderLabel: (key: K, label?: string) => ReactElement;
  showAllNone?: boolean;
}

export const ParamChecklist = <K extends string = string>({
  title,
  checklistState,
  onChecklistState,
  renderLabel,
  options,
  showAllNone = true,
}: PropsWithChildren<ParamChecklistProps<K>>) => {
  const isAll = Object.values(checklistState).every((value) => value);
  const isNone = Object.values(checklistState).every((value) => !value);

  const valueLabels = useMemo(() => options.map(getValueLabel) as ValueLabel<K>[], [options]);
  const keys = valueLabels.map((valueLabel) => valueLabel.value) as K[];

  return (
    <FormGroup>
      <Stack direction="row" alignItems="center" justifyContent="space-between">
        {title && <FormLabel>{title}</FormLabel>}
        {showAllNone && (
          <Box>
            <Button disabled={isAll} onClick={() => onChecklistState(fromKeys(keys, true))}>
              All
            </Button>
            <Button disabled={isNone} onClick={() => onChecklistState(fromKeys(keys, false))}>
              None
            </Button>
          </Box>
        )}
      </Stack>
      {valueLabels.map(({ value: key, label }) => (
        <FormControlLabel
          key={key}
          control={
            <Checkbox
              checked={checklistState[key]}
              onChange={(e, checked) => onChecklistState({ ...checklistState, [key]: checked })}
            />
          }
          label={renderLabel(key, label)}
        />
      ))}
    </FormGroup>
  );
};
