import { ArrowRight } from '@mui/icons-material';
import { Accordion, AccordionDetails, AccordionSummary, Box, Typography } from '@mui/material';
import { FC } from 'react';
import { atomFamily, useRecoilState } from 'recoil';
import { VisibilityToggle } from './VisibilityToggle';

const sidebarSectionExpandedState = atomFamily({
  key: 'sidebarSectionExpandedState',
  default: false,
});

export const SidebarPanel: FC<{ id: string; title: string }> = ({ id, title, children }) => {
  const [expanded, setExpanded] = useRecoilState(sidebarSectionExpandedState(id));

  return (
    <Accordion
      disableGutters
      square // clears the original border radius so that we can set our own
      expanded={expanded}
      onChange={(e, expanded) => setExpanded(expanded)}
      sx={{ pointerEvents: 'auto', marginBottom: 1, borderRadius: 1, overflow: 'hidden' }}
    >
      <AccordionSummary
        sx={{
          '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
            transform: 'rotate(90deg)',
          },
          '& .MuiAccordionSummary-content': {
            marginY: '6px',
          },
          paddingX: '6px',
          flexDirection: 'row-reverse', // this puts the expand icon to the left of the summary bar
        }}
        expandIcon={<ArrowRight />}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', width: '100%' }}>
          <Box sx={{ flexGrow: 1 }}>
            <Typography>{title}</Typography>
          </Box>
          <Box>
            <VisibilityToggle id={id} />
          </Box>
        </Box>
      </AccordionSummary>
      <AccordionDetails sx={{ padding: 0 }}>{children}</AccordionDetails>
    </Accordion>
  );
};
