import { TableCell, TableRow, TextareaAutosize, useTheme } from '@mui/material';
import { AnnotatedValue } from 'config/assessment/effect';
import { Slider, ValueDisplay } from './ValueDisplay';

export const IndicatorRow = ({
  group,
  label,
  strength,
  defaultIndicator,
  revisedIndicator,
  setIndicator,
}: {
  group: string;
  label: string;
  strength: number;
  defaultIndicator: AnnotatedValue;
  revisedIndicator: AnnotatedValue;
  setIndicator: (value: AnnotatedValue) => void;
}) => {
  const theme = useTheme();
  return (
    <>
      <TableRow className={`group-${group}`}>
        <TableCell sx={{ whiteSpace: 'nowrap' }}>{label}</TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <ValueDisplay value={defaultIndicator.value} />
        </TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <Slider
            aria-label={`${label} Revised`}
            value={revisedIndicator.value}
            onChange={(e, value: number) => {
              setIndicator({ ...revisedIndicator, value: value });
            }}
            step={0.01}
            track={false}
            marks={[
              { value: -1, label: '-' },
              { value: 0, label: '0' },
              { value: 1, label: '+' },
            ]}
            min={-1}
            max={1}
            valueLabelDisplay="on"
          />
          <input
            type="number"
            style={{ width: '150px' }}
            value={revisedIndicator.value}
            step={0.01}
            min={-1}
            max={1}
            onChange={(e) => {
              setIndicator({ ...revisedIndicator, value: Number.parseFloat(e.target.value) });
            }}
          />
        </TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <ValueDisplay value={revisedIndicator.value * strength} />
        </TableCell>
      </TableRow>
      {!!revisedIndicator.notes || revisedIndicator.value !== defaultIndicator.value ? (
        <TableRow sx={{ backgroundColor: theme.palette.action.hover }}>
          <TableCell></TableCell>
          <TableCell colSpan={3}>
            <TextareaAutosize
              aria-label="Notes"
              placeholder="Reasons given for change from default effect&hellip;"
              style={{
                width: '100%',
                margin: '5px 0',
                fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
              }}
              minRows={3}
              value={revisedIndicator.notes}
              onChange={(e) => {
                setIndicator({ ...revisedIndicator, notes: e.target.value });
              }}
            />
          </TableCell>
        </TableRow>
      ) : null}
    </>
  );
};
