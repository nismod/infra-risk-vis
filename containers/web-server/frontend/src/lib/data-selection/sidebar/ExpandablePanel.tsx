import { ArrowRight } from '@mui/icons-material';
import {
  AccordionDetailsProps,
  AccordionProps,
  AccordionSummaryProps,
  Box,
  Accordion as MuiAccordion,
  AccordionDetails as MuiAccordionDetails,
  AccordionSummary as MuiAccordionSummary,
  Typography,
} from '@mui/material';
import { styled } from '@mui/styles';
import { FC, ReactChild } from 'react';

export const Accordion = styled(MuiAccordion)({
  pointerEvents: 'auto',
  marginBottom: 1,
  borderRadius: 1,
  overflow: 'hidden',
});

export const AccordionSummary = styled(MuiAccordionSummary)({
  '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
    transform: 'rotate(90deg)',
  },
  flexDirection: 'row-reverse', // this puts the expand icon to the left of the summary bar
  '& .MuiAccordionSummary-content': {
    marginTop: '0',
    marginBottom: '0',
  },
  paddingRight: '5px',
  paddingLeft: '5px',
  paddingTop: '4px',
  paddingBottom: '4px',
  minHeight: '40px',
});

export const AccordionDetails = styled(MuiAccordionDetails)({});

export const AccordionTitle = ({ title, actions }) => {
  return (
    <Box sx={{ display: 'flex', alignItems: 'center', width: '100%' }}>
      <Box sx={{ flexGrow: 1 }}>
        <Typography>{title}</Typography>
      </Box>
      <Box>{actions}</Box>
    </Box>
  );
};

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
      expanded={allowExpand && expanded}
      onChange={(e, expanded) => onExpanded(expanded)}
      disableGutters
      TransitionProps={{ unmountOnExit: true }}
      {...AccordionProps}
    >
      <AccordionSummary
        sx={{
          '& .MuiAccordionSummary-content': {
            marginY: '0',
          },
          paddingX: '0',
        }}
        expandIcon={<ArrowRight color={allowExpand ? 'action' : 'disabled'} />}
        {...AccordionSummaryProps}
      >
        <AccordionTitle title={title} actions={actions} />
      </AccordionSummary>
      <AccordionDetails {...AccordionDetailsProps}>{children}</AccordionDetails>
    </Accordion>
  );
};
