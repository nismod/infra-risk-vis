import { FC } from 'react';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { DroughtsControl } from './DroughtsControl';

export const DroughtsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="drought" title="Drought">
      <SidebarPanelSection>
        <DroughtsControl />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
