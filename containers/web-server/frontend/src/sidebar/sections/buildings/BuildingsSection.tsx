import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { BuildingDensityControl } from './BuildingDensityControl';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <SidebarPanelSection>
        <BuildingDensityControl />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
