import { KeyboardArrowDown, KeyboardArrowUp } from '@mui/icons-material';
import {
  Box,
  Collapse,
  FormControl,
  FormControlLabel,
  FormGroup,
  IconButton,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  Switch,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
} from '@mui/material';
import { AnnotatedValue, Effect } from 'config/assessment/effect';
import { INDICATOR_LABELS } from 'config/assessment/indicators';
import { useState } from 'react';
import { IndicatorRow } from './IndicatorRow';

export function Intervention({
  label,
  defaultEffect,
  revisedEffect,
  options,
  strength,
  setStrength,
  setEffect,
}: {
  label: string;
  defaultEffect: Effect;
  revisedEffect: Effect;
  options: { value: number; label: string }[];
  strength: number;
  setStrength: (e: any) => void;
  setEffect: (key: string, value: AnnotatedValue) => void;
}) {
  const [open, setOpen] = useState(false);
  const [showZeros, setShowZeros] = useState(false);

  const handleZerosSwitch = (e) => {
    setShowZeros(e.target.checked);
  };

  return (
    <>
      {/* Intervention Header Row */}
      <TableRow>
        <TableCell>
          <IconButton aria-label="expand row" size="small" onClick={() => setOpen(!open)}>
            {open ? <KeyboardArrowUp /> : <KeyboardArrowDown />}
          </IconButton>
        </TableCell>
        <TableCell>{label}</TableCell>
        <TableCell>
          <FormControl fullWidth sx={{ my: 1 }}>
            <InputLabel>{label}</InputLabel>
            <Select value={strength} label={label} onChange={setStrength}>
              {options.map(({ value, label }) => (
                <MenuItem key={value} value={value}>
                  {label}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
        </TableCell>
      </TableRow>
      {/* Collapsible Intervention Details */}
      <TableRow>
        <TableCell colSpan={3} sx={{ p: 0 }}>
          <Collapse in={open} timeout="auto" unmountOnExit>
            <Box sx={{ m: 2 }}>
              <FormGroup>
                <FormControlLabel
                  control={<Switch checked={showZeros} onChange={handleZerosSwitch} />}
                  label="Show all indicators"
                />
              </FormGroup>
              <TableContainer component={Paper} sx={{ overflow: 'hidden' }}>
                <Table size="small">
                  <colgroup>
                    <col width="50%" />
                    <col width="15%" />
                    <col width="20%" />
                    <col width="15%" />
                  </colgroup>
                  <TableHead>
                    <TableRow className="group-header">
                      <TableCell>Indicator</TableCell>
                      <TableCell>Default</TableCell>
                      <TableCell>Revised</TableCell>
                      <TableCell>Effective</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {INDICATOR_LABELS.map((option) => {
                      let { value, label } = option;
                      const key = value; // confusingly, take "value" as key into effects objects
                      return revisedEffect &&
                        (showZeros ||
                          defaultEffect[key].value !== 0 ||
                          revisedEffect[key].value !== 0 ||
                          revisedEffect[key].notes) ? (
                        <IndicatorRow
                          key={key}
                          group={key.split('_')[0]}
                          defaultIndicator={defaultEffect[key]}
                          revisedIndicator={revisedEffect[key]}
                          strength={strength}
                          label={label}
                          setIndicator={(indicator) => setEffect(key, indicator)}
                        />
                      ) : null;
                    })}
                  </TableBody>
                </Table>
              </TableContainer>
            </Box>
          </Collapse>
        </TableCell>
      </TableRow>
    </>
  );
}
