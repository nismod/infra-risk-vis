import { ArrowDropUp, ArrowRight } from '@mui/icons-material';
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Box,
  Checkbox,
  FormControlLabel,
  Typography,
} from '@mui/material';
import { FC, useEffect, useState } from 'react';
import { useRecoilState, useSetRecoilState } from 'recoil';
import { showPopulationState } from 'state/population';
import { viewModeState } from 'state/view-mode';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { RegionsControl } from './controls/RegionsControl';
import { ViewModeToggle } from './ViewModeToggle';

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

export const SidebarContent = ({ view }) => {
  const setViewMode = useSetRecoilState(viewModeState);

  useEffect(() => {
    const viewMode = view === 'risk' ? 'direct-damages' : 'input';
    setViewMode(viewMode);
  }, [setViewMode, view]);

  if (view === 'exposure')
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
  else if (view === 'risk')
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
};
