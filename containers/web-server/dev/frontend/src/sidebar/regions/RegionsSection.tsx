import { FC } from 'react';
import { RegionsControl } from './RegionsControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';

export const RegionsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="regions" title="Regions">
      <ErrorBoundary message="There was a problem displaying this section.">
        <SidebarPanelSection>
          <RegionsControl />
        </SidebarPanelSection>
        <SidebarPanelSection variant="style">
          <StyleSelection id="regions" />
        </SidebarPanelSection>
      </ErrorBoundary>
    </SidebarPanel>
  );
};
