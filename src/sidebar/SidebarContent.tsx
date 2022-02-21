import { Box, Typography } from '@mui/material';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { ViewModeToggle } from './ViewModeToggle';

export const SidebarContent = ({ view }) => {
  return (
    view === 'exposure' && (
      <>
        <NetworkControl />
        <Box mb={1}>
          <Typography variant="h6">View Mode</Typography>
          <ViewModeToggle />
        </Box>
        <HazardsControl />
      </>
    )
  );
};
