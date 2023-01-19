import { KeyboardArrowUp, KeyboardArrowDown } from '@mui/icons-material';
import { TableRow, TableCell, IconButton, Collapse, Table, TableBody } from '@mui/material';
import { useState, Fragment } from 'react';
import { useRecoilState } from 'recoil';

import { weightedSum } from 'config/assessment/assessment';
import { Effect } from 'config/assessment/effect';
import { INDICATOR_LABELS } from 'config/assessment/indicators';
import { indicatorWeights } from 'state/assessment';
import { IndicatorTableColGroup } from './IndicatorTableColGroup';
import { ValueDisplay } from './ValueDisplay';
import { WeightDisplay } from './WeightDisplay';
import { WeightRow } from './WeightRow';

export const WeightGroup = ({ label, prefix, unweighted }: { label: string; prefix: string; unweighted: Effect }) => {
  const [open, setOpen] = useState(false);
  const [currentWeights, setWeights] = useRecoilState(indicatorWeights);

  // Force our way around type-checking - recoil returns the expected Effect
  // @ts-ignore
  let weights: Effect = { ...currentWeights };

  const [assessed_value, weighted_value, total_weight] = weightedSum(unweighted, weights, prefix);

  return (
    <Fragment>
      <TableRow className={`group-${prefix}`}>
        <TableCell>
          <IconButton aria-label="expand row" size="small" onClick={() => setOpen(!open)}>
            {open ? <KeyboardArrowUp /> : <KeyboardArrowDown />}
          </IconButton>
        </TableCell>
        <TableCell>{label}</TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <ValueDisplay value={assessed_value} />
        </TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <WeightDisplay value={total_weight} />
        </TableCell>
        <TableCell sx={{ verticalAlign: 'top' }}>
          <ValueDisplay value={weighted_value} />
        </TableCell>
      </TableRow>
      <TableRow className={`group-${prefix}`}>
        <TableCell colSpan={5} sx={{ p: 0 }}>
          <Collapse in={open} timeout="auto" unmountOnExit>
            <Table>
              <IndicatorTableColGroup />
              <TableBody>
                {INDICATOR_LABELS.map((option) => {
                  let { value, label, description } = option;
                  const key = value;
                  if (key.includes(prefix)) {
                    return (
                      <WeightRow
                        key={key}
                        label={label}
                        description={description}
                        assessed_value={unweighted[key]}
                        weight={weights[key]}
                        setWeight={(weight) => {
                          setWeights({
                            ...weights,
                            [key]: weight,
                          });
                        }}
                      />
                    );
                  }
                  return null;
                })}
              </TableBody>
            </Table>
          </Collapse>
        </TableCell>
      </TableRow>
    </Fragment>
  );
};
