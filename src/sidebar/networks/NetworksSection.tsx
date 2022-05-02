import { FC } from 'react';
import { useRecoilValue } from 'recoil';
import { Collapse } from '@mui/material';
import { TransitionGroup } from 'react-transition-group';

import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';

import { NetworkControl } from './NetworkControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { StyleSelection } from 'sidebar/StyleSelection';
import { DamageSourceControl } from './DamageSourceControl';
import { SidebarPanelSection } from 'sidebar/ui/SidebarPanelSection';
import { networksStyleStateEffect, sectionStyleValueState } from 'state/sections';
import { AdaptationControl } from './AdaptationControl';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';

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
