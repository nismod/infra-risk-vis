import { Alert } from '@mui/material';
import { FC, useEffect } from 'react';
import { useSetRecoilState } from 'recoil';

import { viewModeState } from 'state/view-mode';
import { HazardsSection } from './hazards/HazardsSection';
import { NetworksSection } from './networks/NetworksSection';
import { RegionsSection } from './regions/RegionsSection';

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
    case 'risk':
      return (
        <>
          <NetworksSection />
          <HazardsSection />
          <RegionsSection />
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
