import { KeyboardArrowDown, KeyboardArrowUp } from '@mui/icons-material';
import {
  Box,
  Collapse,
  FormControl,
  IconButton,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  Slider,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextareaAutosize,
  useTheme,
} from '@mui/material';
import { AnnotatedValue, Effect } from 'config/assessment/effect';
import { INDICATOR_LABELS } from 'config/assessment/indicators';
import { useState } from 'react';
import { CompactValue } from './CompactValue';
import { ValueDisplay } from './ValueDisplay';

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
  const theme = useTheme();
  return (
    <>
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
            <Slider
              aria-label={label}
              value={strength}
              onChange={(e, value) => {
                setStrength({ target: { value: value } });
              }}
              step={1}
              track={false}
              marks={[
                { value: -1, label: '-' },
                { value: 0, label: '0' },
                { value: 1, label: '+' },
              ]}
              min={-1}
              max={1}
            />
            <CompactValue label="Strength" value={strength} />
          </FormControl>
        </TableCell>
      </TableRow>
      <TableRow>
        <TableCell colSpan={3} sx={{ p: 0 }}>
          <Collapse in={open} timeout="auto" unmountOnExit>
            <Box sx={{ m: 2 }}>
              <TableContainer component={Paper}>
                <Table size="small">
                  <colgroup>
                    <col width="50%" />
                    <col width="15%" />
                    <col width="20%" />
                    <col width="15%" />
                  </colgroup>
                  <TableHead>
                    <TableRow>
                      <TableCell>Indicator</TableCell>
                      <TableCell>Default</TableCell>
                      <TableCell>Revised</TableCell>
                      <TableCell>Effective</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {INDICATOR_LABELS.map((option) => {
                      let { value, label } = option;
                      const key = value;
                      return revisedEffect ? (
                        <>
                          <TableRow key={key}>
                            <TableCell sx={{ whiteSpace: 'nowrap' }}>{label}</TableCell>
                            <TableCell>
                              <ValueDisplay value={defaultEffect[key].value} />
                            </TableCell>
                            <TableCell>
                              <Slider
                                aria-label={`${label} Revised`}
                                value={revisedEffect[key].value}
                                onChange={(e, value: number) => {
                                  setEffect(key, { ...revisedEffect[key], value: value });
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
                              />
                              <input
                                type="number"
                                style={{ width: '150px' }}
                                value={revisedEffect[key].value}
                                step={0.01}
                                min={-1}
                                max={1}
                                onChange={(e) => {
                                  setEffect(key, { ...revisedEffect[key], value: Number.parseFloat(e.target.value) });
                                }}
                              />
                            </TableCell>
                            <TableCell>
                              <ValueDisplay value={revisedEffect[key].value * strength} />
                            </TableCell>
                          </TableRow>
                          {!!revisedEffect[key].notes || revisedEffect[key].value !== defaultEffect[key].value ? (
                            <TableRow key={`annotation-${key}`} sx={{ backgroundColor: theme.palette.action.hover }}>
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
                                  value={revisedEffect[key].notes}
                                  onChange={(e) => {
                                    setEffect(key, { ...revisedEffect[key], notes: e.target.value });
                                  }}
                                />
                              </TableCell>
                            </TableRow>
                          ) : null}
                        </>
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
