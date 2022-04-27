import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { SidebarPanel } from '../SidebarPanel';
import { BuildingsControl } from './BuildingsControl';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <SidebarPanelSection>
        <BuildingsControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection id="buildings" />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
