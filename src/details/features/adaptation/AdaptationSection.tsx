import { Box } from '@mui/system';
import { Typography } from '@mui/material';

import { Adaptation } from 'lib/api-client';
import { AdaptationTable } from './AdaptationTable';


export const AdaptationSection = ({ fd }) => {
  const options : Adaptation[] = fd?.adaptation ?? [];

  return (
    <>
      <Box py={2}>
        <Box position="relative">
          <Typography variant="h6">Adaptation Options</Typography>
          {
            options.length?
              (<AdaptationTable options={options} />)
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
