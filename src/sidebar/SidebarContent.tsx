import { Alert } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { viewState } from 'state/view';

import { HazardsSection } from './hazards/HazardsSection';
import { NetworksSection } from './networks/NetworksSection';
import { RegionsSection } from './regions/RegionsSection';

const viewLabels = {
  exposure: 'Exposure',
  risk: 'Risk',
  adaptation: 'Adaptation',
  'nature-based-solutions': 'Nature-based Solutions',
};

export const SidebarContent: FC<{}> = () => {
  const view = useRecoilValue(viewState);
  switch (view) {
    case 'exposure':
    case 'risk':
    case 'adaptation':
      return (
        <>
          <NetworksSection />
          <HazardsSection />
          {/* <RegionsSection /> */}
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
