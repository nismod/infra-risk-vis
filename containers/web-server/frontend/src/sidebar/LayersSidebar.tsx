import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { SidebarContent } from './SidebarContent';

export const LayersSidebar = () => (
  <ErrorBoundary message="There was a problem displaying the sidebar.">
    <SidebarContent />
  </ErrorBoundary>
);
