import { FC } from 'react';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { MarineControl } from './MarineControl';

export const MarineSection: FC<{}> = () => {
  return (
    <SidebarPanel id="marine" title="Marine">
      <ErrorBoundary message="There was a problem displaying this section.">
        <SidebarPanelSection>
          <MarineControl />
        </SidebarPanelSection>
      </ErrorBoundary>
    </SidebarPanel>
  );
};
