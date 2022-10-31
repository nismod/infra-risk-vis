import { ArrowRight } from '@mui/icons-material';
import { Accordion, AccordionDetails, AccordionSummary, Box, Typography } from '@mui/material';
import { FC, ReactChild } from 'react';

export const ExpandablePanel: FC<{
  expanded: boolean;
  onExpanded: (x: boolean) => void;
  title: string;
  actions: ReactChild;
  disabled?: boolean;
}> = ({ expanded, onExpanded, title, actions, disabled = false, children }) => {
  return (
    <Accordion
      disabled={disabled}
      disableGutters
      square // clears the original border radius so that we can set our own
      expanded={expanded}
      onChange={(e, expanded) => onExpanded(expanded)}
      sx={{
        pointerEvents: 'auto',
        marginBottom: 1,
        borderRadius: 1,
        overflow: 'hidden',
      }}
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
          <Box>{actions}</Box>
        </Box>
      </AccordionSummary>
      <AccordionDetails sx={{ padding: 0 }}>{children}</AccordionDetails>
    </Accordion>
  );
};
