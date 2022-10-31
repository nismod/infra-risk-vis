import { Alert } from '@mui/material';

import { SidebarPanel } from '@/sidebar/SidebarPanel';

export const HealthcareSection = () => {
  return (
    <SidebarPanel id="healthcare" title="Healthcare">
      <Alert severity="info">Use the visibility toggle to show/hide health sites</Alert>
    </SidebarPanel>
  );
};
