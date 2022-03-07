import { FC } from 'react';
import { RegionsControl } from './RegionsControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { regionsStyleState } from 'state/regions';
import { REGION_DEFAULT_STYLE, REGION_STYLES } from 'config/regions/styles';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';

export const RegionsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="regions" title="Regions">
      <SidebarPanelSection>
        <RegionsControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection state={regionsStyleState} options={REGION_STYLES} defaultValue={REGION_DEFAULT_STYLE} />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
