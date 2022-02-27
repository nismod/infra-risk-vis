import { ArrowDropUp, ArrowRight } from '@mui/icons-material';
import { Accordion, AccordionDetails, AccordionSummary, Alert, Typography } from '@mui/material';
import { FC, useEffect, useState } from 'react';
import { useSetRecoilState } from 'recoil';

import { viewModeState } from 'state/view-mode';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { RegionsControl } from './controls/RegionsControl';

const SidebarSection: FC<{ title: string }> = ({ title, children }) => {
  const [expanded, setExpanded] = useState(false);
  return (
    <Accordion expanded={expanded} onChange={(e, expanded) => setExpanded(expanded)} sx={{ pointerEvents: 'auto' }}>
      <AccordionSummary expandIcon={expanded ? <ArrowDropUp /> : <ArrowRight />}>
        <Typography variant="h6">{title}</Typography>
      </AccordionSummary>
      <AccordionDetails>{children}</AccordionDetails>
    </Accordion>
  );
};

const viewLabels = {
  exposure: 'Exposure',
  risk: 'Risk',
  adaptation: 'Adaptation',
  prioritization: 'Prioritization',
};

export const SidebarContent = ({ view }) => {
  const setViewMode = useSetRecoilState(viewModeState);

  useEffect(() => {
    const viewMode = view === 'risk' ? 'direct-damages' : 'input';
    setViewMode(viewMode);
  }, [setViewMode, view]);

  switch (view) {
    case 'exposure':
      return (
        <>
          <SidebarSection title="Built Assets">
            <NetworkControl />
          </SidebarSection>
          <SidebarSection title="Hazards">
            <HazardsControl />
          </SidebarSection>
          <SidebarSection title="Regions">
            <RegionsControl />
          </SidebarSection>
        </>
      );
    case 'risk':
      return (
        <>
          <SidebarSection title="Built Assets">
            <NetworkControl />
          </SidebarSection>
          <SidebarSection title="Hazards">
            <HazardsControl />
          </SidebarSection>
          <SidebarSection title="Regions">
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
