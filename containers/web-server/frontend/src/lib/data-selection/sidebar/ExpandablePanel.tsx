import { ArrowRight } from '@mui/icons-material';
import {
  Accordion,
  AccordionDetails,
  AccordionDetailsProps,
  AccordionProps,
  AccordionSummary,
  AccordionSummaryProps,
  Box,
  Typography,
} from '@mui/material';
import { FC, ReactChild } from 'react';

export const ExpandablePanel: FC<{
  expanded: boolean;
  onExpanded: (x: boolean) => void;
  allowExpand?: boolean;
  title: string;
  actions: ReactChild;
  disabled?: boolean;
  AccordionProps?: Partial<AccordionProps>;
  AccordionSummaryProps?: Partial<AccordionSummaryProps>;
  AccordionDetailsProps?: Partial<AccordionDetailsProps>;
}> = ({
  expanded,
  onExpanded,
  allowExpand = true,
  title,
  actions,
  disabled = false,
  children,
  AccordionProps = {},
  AccordionSummaryProps = {},
  AccordionDetailsProps = {},
}) => {
  return (
    <Accordion
      disabled={disabled}
      disableGutters
      square // clears the original border radius so that we can set our own
      expanded={allowExpand && expanded}
      onChange={(e, expanded) => onExpanded(expanded)}
      sx={{
        pointerEvents: 'auto',
        marginBottom: 1,
        borderRadius: 1,
        overflow: 'hidden',
      }}
      TransitionProps={{ unmountOnExit: true }}
      {...AccordionProps}
    >
      <AccordionSummary
        sx={{
          '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
            transform: 'rotate(90deg)',
          },
          '& .MuiAccordionSummary-content': {
            marginY: '2px',
          },
          paddingX: '2px',
          flexDirection: 'row-reverse', // this puts the expand icon to the left of the summary bar
        }}
        expandIcon={allowExpand ? <ArrowRight /> : null}
        {...AccordionSummaryProps}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', width: '100%' }}>
          <Box sx={{ flexGrow: 1 }}>
            <Typography>{title}</Typography>
          </Box>
          <Box>{actions}</Box>
        </Box>
      </AccordionSummary>
      <AccordionDetails sx={{ padding: 0 }} {...AccordionDetailsProps}>
        {children}
      </AccordionDetails>
    </Accordion>
  );
};
