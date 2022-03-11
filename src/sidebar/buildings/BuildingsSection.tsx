import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { SidebarPanel } from '../SidebarPanel';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <SidebarPanelSection variant="style">
        <StyleSelection id="buildings" />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
