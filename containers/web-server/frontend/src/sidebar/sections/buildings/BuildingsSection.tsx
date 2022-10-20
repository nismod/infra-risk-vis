import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { StyleSelection } from '@/sidebar/StyleSelection';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { BuildingsControl } from './BuildingsControl';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings" disabled>
      <SidebarPanelSection>
        <BuildingsControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection id="buildings" />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
