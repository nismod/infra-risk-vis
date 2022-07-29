import { Box } from '@mui/system';
import { IconButton, Typography } from '@mui/material';

import { Adaptation } from 'lib/api-client';
import { AdaptationTable } from './AdaptationTable';
import { Download } from '@mui/icons-material';
import { downloadFile, unique } from 'lib/helpers';
import { useMemo } from 'react';


function makeAdaptationCsv(options: Adaptation[]) {
  const header = [
    'adaptation_name',
    'hazard',
    'rcp',
    'adaptation_protection_level',
    'adaptation_cost',
    'avoided_ead_amin',
    'avoided_ead_mean',
    'avoided_ead_amax',
    'avoided_eael_amin',
    'avoided_eael_mean',
    'avoided_eael_amax',
  ].join(',') + '\n'

  return header + options.map(
    (d) => [
      `"${d.adaptation_name}"`,
      d.hazard,
      d.rcp,
      d.adaptation_protection_level,
      d.adaptation_cost,
      d.avoided_ead_amin,
      d.avoided_ead_mean,
      d.avoided_ead_amax,
      d.avoided_eael_amin,
      d.avoided_eael_mean,
      d.avoided_eael_amax,
    ].join(',')
  ).join('\n');
}

export const AdaptationSection = ({ fd }) => {
  const options : Adaptation[] = useMemo(() => {
    // TODO: remove factor when data is updated
    const HORRIBLE_HACK_FACTOR = 1/15;
    return fd?.adaptation.map(
      (d: Adaptation) => {
        return {
          adaptation_name: d.adaptation_name,
          hazard: d.hazard,
          rcp: d.rcp,
          adaptation_protection_level: d.adaptation_protection_level,
          adaptation_cost: d.adaptation_cost,
          avoided_ead_amin: d.avoided_ead_amin,
          avoided_ead_mean: d.avoided_ead_mean,
          avoided_ead_amax: d.avoided_ead_amax,
          avoided_eael_amin: d.avoided_eael_amin * HORRIBLE_HACK_FACTOR,
          avoided_eael_mean: d.avoided_eael_mean * HORRIBLE_HACK_FACTOR,
          avoided_eael_amax: d.avoided_eael_amax * HORRIBLE_HACK_FACTOR,
        }
      }
    ) ?? [];
  }, [fd])
  const option_names = useMemo(() => unique(
    options.map((d) => d.adaptation_name)), [options]
  ).sort();

  return (
    <>
      <Box py={2}>
        <Box position="relative">
          <Typography variant="h6">Adaptation Options</Typography>
          <IconButton
            sx={{
              position: 'absolute',
              top: 0,
              right: 0,
            }}
            title="Download CSV with adaptation options data"
            onClick={() => downloadFile(makeAdaptationCsv(options), 'text/csv', `feature_${fd.id}_adaptation.csv`)}
          >
            <Download />
          </IconButton>
        </Box>
        <Box>
          {
            options.length?
              option_names.map(name => {
                const filteredOptions = options.filter(o=>o.adaptation_name === name);
                return (
                  <>
                    <Typography variant="subtitle2" component="h2" sx={{mt:2,mb:1}}>{name}</Typography>
                    <Typography variant="body2" component="p" sx={{mb:2}}>
                      The adaptation costs and benefits, subject to different
                      climate scenarios and (for some options) protection
                      standards, assuming a 15-day disruption.
                    </Typography>

                    <AdaptationTable options={filteredOptions} />
                  </>
                )
              })
              :
              (
                <Typography variant="body2" color="textSecondary">
                  No adaptation options evaluated.
                </Typography>
              )
          }
        </Box>
      </Box>
    </>
  )
}
