import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';

import { numFormat, numRangeFormat } from 'lib/helpers';

const padding =  {px:0.25,py:0.25}

export const AdaptationTable = ({ options }) => {

  const numDays = 15;
  return (
    <TableContainer sx={{maxHeight:260}}>
      <Table size="small" padding="none" stickyHeader>
        <TableHead>
          <TableRow>
            <TableCell sx={padding}>Hazard</TableCell>
            <TableCell sx={padding}>
              <abbr title="Representative Concentration Pathway (Climate Scenario)">RCP</abbr>
            </TableCell>
            <TableCell sx={padding}>Protection Standard</TableCell>
            <TableCell sx={padding} align="right">
              <abbr title="Benefit Cost Ratio">BCR</abbr>
            </TableCell>
            <TableCell sx={padding} align="right">Cost</TableCell>
            <TableCell sx={{pr:0,pl:padding.px,py:padding.py}} align="right">
              Avoided Risk (J$)
            </TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {options.map((d, i) => (
            <TableRow>
              <TableCell sx={{pl:0,pr:padding.px,py:padding.py}}>{d.hazard}</TableCell>
              <TableCell sx={padding}>{d.rcp}</TableCell>
              <TableCell sx={padding}>{d.adaptation_protection_level}</TableCell>
              <TableCell sx={padding} align="right">
                {
                  numFormat(
                    (d.avoided_ead_mean + (d.avoided_eael_mean * numDays)) / d.adaptation_cost
                  )
                }<br/>{
                  numRangeFormat(
                    (d.avoided_ead_amin + (d.avoided_eael_amin * numDays)) / d.adaptation_cost,
                    (d.avoided_ead_amax + (d.avoided_eael_amax * numDays)) / d.adaptation_cost
                  )
                }
              </TableCell>
              <TableCell sx={padding} align="right">{numFormat(d.adaptation_cost)}</TableCell>
              <TableCell sx={{pr:0,pl:padding.px,py:padding.py}} align="right">
                {
                  numFormat(d.avoided_ead_mean + (d.avoided_eael_mean * numDays))
                }<br/>{
                  numRangeFormat(
                    (d.avoided_ead_amin + (d.avoided_eael_amin * numDays)),
                    (d.avoided_ead_amax + (d.avoided_eael_amax * numDays))
                  )
                }
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  )
}
