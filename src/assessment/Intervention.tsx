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
  Slider,
  Switch,
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

const IndicatorRow = ({
  label,
  strength,
  defaultIndicator,
  revisedIndicator,
  setIndicator,
}: {
  label: string;
  strength: number;
  defaultIndicator: AnnotatedValue;
  revisedIndicator: AnnotatedValue;
  setIndicator: (value: AnnotatedValue) => void;
}) => {
    const theme = useTheme();
    return (
      <>
        <TableRow>
          <TableCell sx={{ whiteSpace: 'nowrap' }}>{label}</TableCell>
          <TableCell>
            <ValueDisplay value={defaultIndicator.value} />
          </TableCell>
          <TableCell>
            <Slider
              aria-label={`${label} Revised`}
              value={revisedIndicator.value}
              onChange={(e, value: number) => {
                setIndicator({ ...revisedIndicator, value: value });
              } }
              step={0.01}
              track={false}
              marks={[
                { value: -1, label: '-' },
                { value: 0, label: '0' },
                { value: 1, label: '+' },
              ]}
              min={-1}
              max={1} />
            <input
              type="number"
              style={{ width: '150px' }}
              value={revisedIndicator.value}
              step={0.01}
              min={-1}
              max={1}
              onChange={(e) => {
                setIndicator({ ...revisedIndicator, value: Number.parseFloat(e.target.value) });
              } } />
          </TableCell>
          <TableCell>
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
                } } />
            </TableCell>
          </TableRow>
        ) : null}
      </>
    );
  };

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
      {/* Collapsible Intervention Details */}
      <TableRow>
        <TableCell colSpan={3} sx={{ p: 0 }}>
          <Collapse in={open} timeout="auto" unmountOnExit>
            <Box sx={{ m: 2 }}>
              <FormGroup>
                <FormControlLabel control={<Switch checked={showZeros} onChange={handleZerosSwitch} />} label="Show zero-value indicators" />
              </FormGroup>
              <TableContainer component={Paper}>
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
                      const key = value;  // confusingly, take "value" as key into effects objects
                      return (revisedEffect && (showZeros || (defaultEffect[key].value !== 0 || revisedEffect[key].value !== 0 || revisedEffect[key].notes))) ? 
                        <IndicatorRow 
                          key={key} 
                          group={key.split("_")[0]}
                          defaultIndicator={defaultEffect[key]}
                          revisedIndicator={revisedEffect[key]}
                          strength={strength}
                          label={label}
                          setIndicator={(indicator)=> setEffect(key, indicator)}
                        /> 
                        : null;
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
