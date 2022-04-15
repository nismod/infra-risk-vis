import { Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import { numFormat } from 'lib/helpers';

export const DamageTable = ({ damages }) => (
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
        {damages.map(({ key, hazard, rcp, epoch, ead }) => (
          <TableRow key={key}>
            <TableCell>{hazard}</TableCell>
            <TableCell>{rcp}</TableCell>
            <TableCell>{epoch}</TableCell>
            <TableCell align="right">{numFormat(ead)}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  </TableContainer>
);
