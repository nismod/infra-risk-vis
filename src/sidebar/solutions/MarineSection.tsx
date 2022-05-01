import { FC } from 'react';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { MarineControl } from './MarineControl';

export const MarineSection: FC<{}> = () => {
  return (
    <SidebarPanel id="marine" title="Marine Solutions">
      <SidebarPanelSection>
        <MarineControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection id="marine" />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
