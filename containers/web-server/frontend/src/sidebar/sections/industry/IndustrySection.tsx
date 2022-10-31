import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { IndustryControl } from './IndustryControl';

export const IndustrySection = () => {
  return (
    <SidebarPanel id="industry" title="Industry">
      <SidebarPanelSection>
        <IndustryControl />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
