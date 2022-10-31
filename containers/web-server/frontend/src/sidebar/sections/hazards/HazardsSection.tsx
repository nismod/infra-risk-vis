import { Box } from '@mui/system';
import { FC } from 'react';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { SidebarPanel } from '@/sidebar/SidebarPanel';

import { HazardsControl } from './HazardsControl';

export const HazardsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="hazards" title="Hazards">
      <Box p={2}>
        <ErrorBoundary message="There was a problem displaying this section.">
          <HazardsControl />
        </ErrorBoundary>
      </Box>
    </SidebarPanel>
  );
};
