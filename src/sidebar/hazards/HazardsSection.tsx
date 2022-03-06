import { FC } from 'react';
import { HazardsControl } from './HazardsControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { Box } from '@mui/system';

export const HazardsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="hazards" title="Hazards">
      <Box p={2}>
        <HazardsControl />
      </Box>
    </SidebarPanel>
  );
};
