import React, { FC } from 'react';
import Typography from '@material-ui/core/Typography';
import { Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@material-ui/core';


import { hazardConfig } from '../controls/use-hazard-selection';
import { numFormat } from '../helpers';

interface RiskSectionProps {
  f: any
}

export const RiskSection: FC<RiskSectionProps> = ({ f }) => {
  const rows = Object.entries(hazardConfig).map(([hazardType, hazard]) => {
    const items = [];
    for (const rcp of hazard.paramDomains.rcp) {
      for (const epoch of hazard.paramDomains.epoch) {
        // TODO check risk data for confidence
        // for (const confidence of hazard.paramDomains.confidence) {
        const risk_key = `${hazardType}__rcp_${rcp}__epoch_${epoch}__conf_None`
        if (f[risk_key]) {
          items.push((
          <TableRow key={risk_key}>
            <TableCell>{hazardType}</TableCell>
            <TableCell>{rcp}</TableCell>
            <TableCell>{epoch}</TableCell>
            <TableCell align="right">{numFormat(f[risk_key])}</TableCell>
          </TableRow>
          ))
        }
        // }
      }
    }
    return items;
  }).flat();

  return (
    <>
    <Typography variant="subtitle2">Risk</Typography>
    {
      rows.length?
      (
        <TableContainer component={Paper}>
          <Table>
            <TableHead>
              <TableRow>
                <TableCell>Hazard</TableCell>
                <TableCell><abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr></TableCell>
                <TableCell>Epoch</TableCell>
                <TableCell align="right"><abbr title="Expected Annual Damages">EAD</abbr></TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
            {
              rows
            }
            </TableBody>
          </Table>
        </TableContainer>
      )
      : <Typography
          variant="body2"
          color="textSecondary"
        >No exposure direct damages estimated.</Typography>
    }
    </>
  );
}
