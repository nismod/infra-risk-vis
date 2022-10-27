import { FC } from 'react';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { StyleSelection } from '@/sidebar/StyleSelection';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';

import { TerrestrialControl } from './TerrestrialControl';

export const TerrestrialSection: FC<{}> = () => {
  return (
    <SidebarPanel id="terrestrial" title="Terrestrial">
      <SidebarPanelSection>
        <TerrestrialControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection id="terrestrial" />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
