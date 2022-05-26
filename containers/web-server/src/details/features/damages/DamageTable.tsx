import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import { numFormat, numRangeFormat } from 'lib/helpers';
const padding =  {px:0.25,py:0.25}
export const DamageTable = ({ damages }) => (
  <TableContainer sx={{maxHeight:260}}>
    <Table size="small" padding="none" stickyHeader>
      <TableHead>
        <TableRow>
          <TableCell sx={padding}>
            <abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr>
          </TableCell>
          <TableCell sx={padding}>Epoch</TableCell>
          <TableCell sx={padding} align="right">
            <abbr title="Expected Annual Damages">EAD (J$)</abbr>
          </TableCell>
          <TableCell sx={{pr:0,pl:padding.px,py:padding.py}} align="right">
            <abbr title="Expected Annual Economic Losses">EAEL</abbr> (J$/day)
          </TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {damages.map(({ key, rcp, epoch, ead_mean, ead_amin, ead_amax, eael_mean, eael_amin, eael_amax }) => (
          <TableRow key={key}>
            <TableCell sx={{pl:0,pr:padding.px,py:padding.py}}>{rcp}</TableCell>
            <TableCell sx={padding}>{epoch}</TableCell>
            <TableCell sx={padding} align="right">
              {numFormat(ead_mean)}<br/>({numRangeFormat(ead_amin, ead_amax)})
            </TableCell>
            <TableCell sx={{pr:0,pl:padding.px,py:padding.py}} align="right">
              {
                eael_amax? numFormat(eael_mean) : '-'
              }<br/>{
                eael_amax? `(${numRangeFormat(eael_amin, eael_amax)})` : null
              }

            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  </TableContainer>
);
