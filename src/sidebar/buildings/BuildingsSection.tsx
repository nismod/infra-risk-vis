import { StyleSelection } from 'sidebar/StyleSelection';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { buildingsStyleState } from 'state/buildings';
import { SidebarPanel } from '../SidebarPanel';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <SidebarPanelSection variant="style">
        <StyleSelection
          state={buildingsStyleState}
          options={[
            {
              label: 'Building type',
              value: 'type',
            },
          ]}
          defaultValue="type"
        />
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
