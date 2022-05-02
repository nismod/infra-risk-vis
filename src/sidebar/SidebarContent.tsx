import { Alert } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { viewState } from 'state/view';

import { BuildingsSection } from './buildings/BuildingsSection';
import { DroughtsSection } from './drought/DroughtsSection';
import { HazardsSection } from './hazards/HazardsSection';
import { NetworksSection } from './networks/NetworksSection';
import { RegionsSection } from './regions/RegionsSection';
import { MarineSection } from './solutions/MarineSection';
import { TerrestrialSection } from './solutions/TerrestrialSection';

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
      return (
        <>
          <NetworksSection />
          <HazardsSection />
          <BuildingsSection />
          <RegionsSection />
        </>
      );
    case 'adaptation':
      return (
        <>
          <NetworksSection />
          <DroughtsSection />
          <HazardsSection />
          <BuildingsSection />
          <RegionsSection />
          <TerrestrialSection />
          <MarineSection />
        </>
      );
    case 'nature-based-solutions':
      return (
        <>
          <TerrestrialSection />
          <MarineSection />
          <NetworksSection />
          <HazardsSection />
          <BuildingsSection />
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
