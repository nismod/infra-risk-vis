import { FC } from 'react';
import { NetworkControl } from './NetworkControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { DamageSourceControl } from './DamageSourceControl';
import { useRecoilValue } from 'recoil';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { sectionStyleValueState } from 'state/sections';
import { TransitionGroup } from 'react-transition-group';
import { Collapse } from '@mui/material';

export const NetworksSection: FC<{}> = () => {
  const style = useRecoilValue(sectionStyleValueState('assets'));
  return (
    <SidebarPanel id="assets" title="Infrastructure">
      <SidebarPanelSection>
        <NetworkControl />
      </SidebarPanelSection>
      <SidebarPanelSection variant="style">
        <StyleSelection id="assets" />
        <TransitionGroup>
          {style === 'direct-damages' ? (
            <Collapse>
              <DamageSourceControl />
            </Collapse>
          ) : null}
        </TransitionGroup>
      </SidebarPanelSection>
    </SidebarPanel>
  );
};
