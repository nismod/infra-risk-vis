import { Alert } from '@mui/material';
import { FC, useEffect } from 'react';
import { useSetRecoilState } from 'recoil';

import { viewModeState } from 'state/view-mode';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { RegionsControl } from './controls/RegionsControl';
import { SidebarSection } from './SidebarSection';

const viewLabels = {
  exposure: 'Exposure',
  risk: 'Risk',
  adaptation: 'Adaptation',
  prioritization: 'Prioritization',
};

export const SidebarContent: FC<{ view: string }> = ({ view }) => {
  const setViewMode = useSetRecoilState(viewModeState);

  useEffect(() => {
    const viewMode = view === 'risk' ? 'direct-damages' : 'input';
    setViewMode(viewMode);
  }, [setViewMode, view]);

  switch (view) {
    case 'exposure':
      return (
        <>
          <SidebarSection id="assets" title="Built Assets">
            <NetworkControl />
          </SidebarSection>
          <SidebarSection id="hazards" title="Hazards">
            <HazardsControl />
          </SidebarSection>
          <SidebarSection id="regions" title="Regions">
            <RegionsControl />
          </SidebarSection>
        </>
      );
    case 'risk':
      return (
        <>
          <SidebarSection id="assets" title="Built Assets">
            <NetworkControl />
          </SidebarSection>
          <SidebarSection id="hazards" title="Hazards">
            <HazardsControl />
          </SidebarSection>
          <SidebarSection id="regions" title="Regions">
            <RegionsControl />
          </SidebarSection>
        </>
      );
    default: {
      const viewLabel = viewLabels[view];

      if (viewLabel) {
        return <Alert severity="info">{viewLabel}: Coming soon.</Alert>;
      } else {
        return <Alert severity="error">Unknown view!</Alert>;
      }
    }
  }
};
