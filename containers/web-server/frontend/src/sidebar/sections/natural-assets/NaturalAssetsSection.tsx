import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { NaturalAssetsControl } from './NaturalAssetsControl';

export const NaturalAssetsSection = () => {
  return (
    <SidebarPanel id="natural-assets" title="Natural Assets">
      <SidebarPanelSection>
        <NaturalAssetsControl />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
