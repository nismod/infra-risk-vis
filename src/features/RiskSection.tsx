import React, { FC } from 'react';
import Typography from '@material-ui/core/Typography';
import { Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@material-ui/core';

import { numFormat } from '../helpers';

interface RiskSectionProps {
  eadData: any;
}

export const RiskSection: FC<RiskSectionProps> = ({ eadData }) => {
  return (
    <>
      <Typography variant="subtitle2">Risk</Typography>
      {eadData.length ? (
        <TableContainer component={Paper}>
          <Table>
            <TableHead>
              <TableRow>
                <TableCell>Hazard</TableCell>
                <TableCell>
                  <abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr>
                </TableCell>
                <TableCell>Epoch</TableCell>
                <TableCell align="right">
                  <abbr title="Expected Annual Damages">EAD</abbr>
                </TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {eadData.map(({ key, hazardType, rcp, epoch, ead }) => (
                <TableRow key={key}>
                  <TableCell>{hazardType}</TableCell>
                  <TableCell>{rcp}</TableCell>
                  <TableCell>{epoch}</TableCell>
                  <TableCell align="right">{numFormat(ead)}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      ) : (
        <Typography variant="body2" color="textSecondary">
          No exposure direct damages estimated.
        </Typography>
      )}
    </>
  );
};
