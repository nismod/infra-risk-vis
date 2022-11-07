import { Box } from '@mui/system';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { DeselectButton } from './DeselectButton';
import { SidePanel } from './SidePanel';

export const DetailsPanel = ({ interactionGroup, children }) => {
  return (
    <SidePanel position="relative">
      <ErrorBoundary message="There was a problem displaying these details.">
        <Box position="absolute" top={0} right={0} p={2} zIndex={1000}>
          <DeselectButton interactionGroup={interactionGroup} title="Deselect" />
        </Box>
        {children}
      </ErrorBoundary>
    </SidePanel>
  );
};
