import { Alert } from '@mui/material';

import { SidebarPanel } from '@/sidebar/SidebarPanel';

export const PopulationSection = () => {
  return (
    <SidebarPanel id="population" title="Population">
      <Alert severity="info">Use the visibility toggle to show/hide population</Alert>
    </SidebarPanel>
  );
};
