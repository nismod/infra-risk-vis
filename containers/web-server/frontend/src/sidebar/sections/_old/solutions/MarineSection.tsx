import { FC } from 'react';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { MarineControl } from './MarineControl';

export const MarineSection: FC<{}> = () => {
  return (
    <SidebarPanel id="marine" title="Marine">
      <SidebarPanelSection>
        <MarineControl />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
