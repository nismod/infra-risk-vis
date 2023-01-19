import { TableRow, TableCell } from '@mui/material';
import { AnnotatedValue } from 'config/assessment/effect';
import { HelpTooltip } from './HelpTooltip';
import { Slider, ValueDisplay } from './ValueDisplay';

export const WeightRow = ({
  label,
  description,
  assessed_value,
  weight,
  setWeight,
}: {
  label: string;
  description: string;
  assessed_value: AnnotatedValue;
  weight: AnnotatedValue;
  setWeight: (value: AnnotatedValue) => void;
}) => {
  return (
    <TableRow>
      <TableCell />
      <TableCell sx={{ verticalAlign: 'top', py: 2 }}>
        {label}
        {HelpTooltip(description)}
      </TableCell>
      <TableCell sx={{ verticalAlign: 'top' }}>
        <ValueDisplay value={assessed_value.value} />
      </TableCell>
      <TableCell sx={{ verticalAlign: 'top' }}>
        <Slider
          aria-label={label}
          value={weight.value}
          onChange={(e, value: number) => {
            setWeight({ ...weight, value: value });
          }}
          step={0.01}
          marks={[
            { value: 0, label: '0' },
            { value: 1, label: '+' },
          ]}
          min={0}
          max={1}
          valueLabelDisplay="on"
        />
        <input
          type="number"
          style={{ width: '150px' }}
          value={weight.value}
          step={0.01}
          min={0}
          max={1}
          onChange={(e) => {
            setWeight({ ...weight, value: Number.parseFloat(e.target.value) });
          }}
        />
      </TableCell>
      <TableCell sx={{ verticalAlign: 'top' }}>
        <ValueDisplay value={assessed_value.value * weight.value} />
      </TableCell>
    </TableRow>
  );
};
