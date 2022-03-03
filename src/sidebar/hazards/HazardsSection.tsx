import { FC } from 'react';
import { HazardsControl } from './HazardsControl';
import { SidebarSection } from 'sidebar/SidebarSection';

export const HazardsSection: FC<{}> = () => {
  return (
    <SidebarSection id="hazards" title="Hazards">
      <HazardsControl />
    </SidebarSection>
  );
};
