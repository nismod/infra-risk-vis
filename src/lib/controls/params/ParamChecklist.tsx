import { Button, Checkbox, FormControlLabel, FormGroup, FormLabel, Stack } from '@mui/material';
import { Box } from '@mui/system';
import { fromKeys } from 'lib/helpers';
import { PropsWithChildren, ReactElement } from 'react';

interface ParamChecklistProps<K extends string> {
  title: string;
  options: K[];
  checklistState: { [key in K]: boolean };
  onChecklistState: (state: { [key in K]: boolean }) => void;
  renderLabel: (key: K) => ReactElement;
}

export const ParamChecklist = <K extends string = string>({
  title,
  checklistState,
  onChecklistState,
  renderLabel,
  options,
}: PropsWithChildren<ParamChecklistProps<K>>) => {
  const isAll = Object.values(checklistState).every((value) => value);
  const isNone = Object.values(checklistState).every((value) => !value);
  return (
    <FormGroup>
      <Stack direction="row" alignItems="center" justifyContent="space-between">
        <FormLabel>{title}</FormLabel>
        <Box>
          <Button disabled={isAll} onClick={() => onChecklistState(fromKeys(options, true))}>
            All
          </Button>
          <Button disabled={isNone} onClick={() => onChecklistState(fromKeys(options, false))}>
            None
          </Button>
        </Box>
      </Stack>
      {options
        .map((o) => [o, checklistState[o]] as [K, boolean])
        .map(([key, value]) => (
          <FormControlLabel
            key={key}
            control={
              <Checkbox
                checked={value}
                onChange={(e, checked) => onChecklistState({ ...checklistState, [key]: checked })}
              />
            }
            label={renderLabel(key)}
          />
        ))}
    </FormGroup>
  );
};
