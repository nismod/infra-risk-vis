import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { SidebarPanel } from '../SidebarPanel';
import { BuildingsControl } from './BuildingsControl';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <ErrorBoundary message="There was a problem displaying this section.">
        <SidebarPanelSection>
          <BuildingsControl />
        </SidebarPanelSection>
        <SidebarPanelSection variant="style">
          <StyleSelection id="buildings" />
        </SidebarPanelSection>
      </ErrorBoundary>
    </SidebarPanel>
  );
};
