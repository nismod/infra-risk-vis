import { FC } from 'react';

import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';
import { networksStyleStateEffect, sectionStyleValueState } from '@/state/sections';

import { NetworkControl } from './NetworkControl';

export const NetworksSection: FC<{}> = () => {
  // const style = useRecoilValue(sectionStyleValueState('assets'));
  return (
    <SidebarPanel id="assets" title="Infrastructure">
      <StateEffectRoot state={sectionStyleValueState('assets')} effect={networksStyleStateEffect} />
      <SidebarPanelSection>
        <NetworkControl />
      </SidebarPanelSection>
      {/* <SidebarPanelSection variant="style">
        <StyleSelection id="assets" />
        <TransitionGroup>
          {style === 'damages' ? (
            <Collapse>
              <DamageSourceControl />
            </Collapse>
          ) : null}
          {style === 'adaptation' ? (
            <Collapse>
              <AdaptationControl />
            </Collapse>
          ) : null}
        </TransitionGroup>
      </SidebarPanelSection> */}
    </SidebarPanel>
  );
};
