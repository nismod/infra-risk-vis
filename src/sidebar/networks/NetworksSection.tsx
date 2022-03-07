import { FC } from 'react';
import { NetworkControl } from './NetworkControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { NETWORKS_DEFAULT_STYLE, NETWORK_STYLES } from 'config/networks/styles';
import { networksStyleState } from 'state/networks/networks-style';
import { DamageSourceControl } from './DamageSourceControl';
import { useRecoilValue } from 'recoil';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';

export const NetworksSection: FC<{}> = () => {
  const style = useRecoilValue(networksStyleState);
  return (
    <SidebarPanel id="assets" title="Infrastructure">
      <SidebarPanelSection>
        <NetworkControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection state={networksStyleState} options={NETWORK_STYLES} defaultValue={NETWORKS_DEFAULT_STYLE} />
        {style === 'direct-damages' ? <DamageSourceControl /> : null}
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
