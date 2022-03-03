import { FC } from 'react';
import { RegionsControl } from './RegionsControl';
import { SidebarSection } from 'sidebar/SidebarSection';

export const RegionsSection: FC<{}> = () => {
  return (
    <SidebarSection id="regions" title="Regions">
      <RegionsControl />
    </SidebarSection>
  );
};
