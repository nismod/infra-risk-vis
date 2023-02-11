import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { MobileTabContentWatcher } from '@/pages/map/layouts/mobile/tab-has-content';

import { SidebarContent } from './SidebarContent';

export const LayersSidebar = () => (
  <ErrorBoundary message="There was a problem displaying the sidebar.">
    <MobileTabContentWatcher tabId="layers" />
    <SidebarContent />
  </ErrorBoundary>
);
