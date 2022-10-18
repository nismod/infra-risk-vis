import { Collapse } from '@mui/material';
import { FC } from 'react';
import { TransitionGroup } from 'react-transition-group';
import { useRecoilValue } from 'recoil';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { SidebarPanel } from '@/sidebar/SidebarPanel';
import { StyleSelection } from '@/sidebar/StyleSelection';
import { SidebarPanelSection } from '@/sidebar/ui/SidebarPanelSection';
import { networksStyleStateEffect, sectionStyleValueState } from '@/state/sections';

import { AdaptationControl } from './AdaptationControl';
import { DamageSourceControl } from './DamageSourceControl';
import { NetworkControl } from './NetworkControl';

export const NetworksSection: FC<{}> = () => {
  const style = useRecoilValue(sectionStyleValueState('assets'));
  return (
    <SidebarPanel id="assets" title="Infrastructure">
      <ErrorBoundary message="There was a problem displaying this section.">
        <StateEffectRoot state={sectionStyleValueState('assets')} effect={networksStyleStateEffect} />
        <SidebarPanelSection>
          <NetworkControl />
        </SidebarPanelSection>
        <SidebarPanelSection variant="style">
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
        </SidebarPanelSection>
      </ErrorBoundary>
    </SidebarPanel>
  );
};
