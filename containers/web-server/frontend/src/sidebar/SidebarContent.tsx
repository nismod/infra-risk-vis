import { Alert } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { viewState } from '@/state/view';

import { BuildingsSection } from './sections/buildings/BuildingsSection';
import { HazardsSection } from './sections/hazards/HazardsSection';
import { HealthcareSection } from './sections/healthcare/HealthcareSection';
import { IndustrySection } from './sections/industry/IndustrySection';
import { NaturalAssetsSection } from './sections/natural-assets/NaturalAssetsSection';
import { NetworksSection } from './sections/networks/NetworksSection';
import { PopulationSection } from './sections/population/PopulationSection';
import { HumanVulnerabilitySection } from './sections/vulnerability/HumanVulnerabilitySection';
import { NatureVulnerabilitySection } from './sections/vulnerability/NatureVulnerabilitySection';

const viewLabels = {
  hazard: 'Hazard',
  exposure: 'Exposure',
  vulnerability: 'Vulnerability',
  risk: 'Risk',
};

export const SidebarContent: FC<{}> = () => {
  const view = useRecoilValue(viewState);
  switch (view) {
    case 'hazard':
    case 'exposure':
    case 'vulnerability':
    case 'risk':
      return (
        <>
          {/* Hazard sections */}
          <HazardsSection />

          {/* Exposure sections */}
          <PopulationSection />
          <BuildingsSection />
          <NetworksSection />
          <IndustrySection />
          <HealthcareSection />
          <NaturalAssetsSection />

          {/* Vulnerability sections */}
          <HumanVulnerabilitySection />
          <NatureVulnerabilitySection />
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
