import { ArrowRight, Visibility, VisibilityOff } from '@mui/icons-material';
import { Accordion, AccordionDetails, AccordionSummary, Alert, IconButton, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { FC, useEffect } from 'react';
import { atomFamily, useRecoilState, useSetRecoilState } from 'recoil';
import { showLayerState } from 'state/show-layers';

import { viewModeState } from 'state/view-mode';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { RegionsControl } from './controls/RegionsControl';

const sidebarSectionExpandedState = atomFamily({
  key: 'sidebarSectionExpandedState',
  default: false,
});

const SidebarSection: FC<{ id: string; title: string }> = ({ id, title, children }) => {
  const [expanded, setExpanded] = useRecoilState(sidebarSectionExpandedState(id));
  const [visibility, setVisibility] = useRecoilState(showLayerState(id));

  return (
    <Accordion
      disableGutters
      square // clears the original border radius so that we can set our own
      expanded={expanded}
      onChange={(e, expanded) => setExpanded(expanded)}
      sx={{ pointerEvents: 'auto', marginBottom: 1, borderRadius: 1 }}
    >
      <AccordionSummary
        sx={{
          '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
            transform: 'rotate(90deg)',
          },
          flexDirection: 'row-reverse', // this puts the expand icon to the left of the summary bar
        }}
        expandIcon={<ArrowRight />}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', width: '100%' }}>
          <Box sx={{ flexGrow: 1 }}>
            <Typography variant="h6">{title}</Typography>
          </Box>
          <Box>
            <IconButton
              onClick={(e) => {
                setVisibility((visibility) => !visibility);
                e.stopPropagation();
              }}
            >
              {visibility ? <Visibility /> : <VisibilityOff />}
            </IconButton>
          </Box>
        </Box>
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
