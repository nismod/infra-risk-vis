import { FC } from 'react';
import { NetworkControl } from './NetworkControl';
import { SidebarSection } from 'sidebar/SidebarSection';

export const NetworksSection: FC<{}> = () => {
  return (
    <SidebarSection id="assets" title="Infrastructure">
      <NetworkControl />
    </SidebarSection>
  );
};
