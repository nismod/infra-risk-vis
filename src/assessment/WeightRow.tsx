import { TableRow, TableCell, Slider } from "@mui/material";
import { ValueDisplay } from "./ValueDisplay";

export const WeightRow = ({ label, assessed_value, weight, setWeight }) => {
  return (
    <TableRow>
      <TableCell />
      <TableCell>{label}</TableCell>
      <TableCell>
        <ValueDisplay value={assessed_value} />
      </TableCell>
      <TableCell>
        <Slider
          aria-label={label}
          value={weight}
          onChange={setWeight}
          step={0.01}
          marks={[
            { value: 0, label: '0' },
            { value: 1, label: '+' },
          ]}
          min={0}
          max={1}
        />
        <input
          type="number"
          style={{ width: '150px' }}
          value={weight}
          step={0.01}
          min={0}
          max={1}
          onChange={(e) => {
            setWeight(e, Number.parseFloat(e.target.value));
          }}
        />
      </TableCell>
      <TableCell>
        <ValueDisplay value={assessed_value * weight} />
      </TableCell>
    </TableRow>
  );
};