import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';

import { numFormat, numRangeFormat } from '@/lib/helpers';

const padding = { px: 0.25, py: 0.25 };
export const RPDamageTable = ({ damages }) => (
  <TableContainer sx={{ maxHeight: 250 }}>
    <Table size="small" padding="none" stickyHeader>
      <TableHead>
        <TableRow>
          <TableCell sx={padding}>
            <abbr title="Return Period">RP</abbr>
          </TableCell>
          <TableCell sx={padding}>
            <abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr>
          </TableCell>
          <TableCell sx={padding} align="right">
            Damages (USD)
          </TableCell>
          <TableCell sx={{ pr: 0, pl: padding.px, py: padding.py }} align="right">
            Loss (USD/day)
          </TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {damages.map((d) => (
          <TableRow key={d.key}>
            <TableCell sx={padding}>{d.rp}</TableCell>
            <TableCell sx={{ pl: 0, pr: padding.px, py: padding.py }}>{d.rcp}</TableCell>
            <TableCell sx={padding} align="right">
              {numFormat(d.damage_mean)}
              <br />({numRangeFormat(d.damage_amin, d.damage_amax)})
            </TableCell>
            <TableCell sx={{ pr: 0, pl: padding.px, py: padding.py }} align="right">
              {d.loss_amax ? numFormat(d.loss_mean) : '-'}
              <br />
              {d.loss_amax ? `(${numRangeFormat(d.loss_amin, d.loss_amax)})` : null}
            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  </TableContainer>
);
