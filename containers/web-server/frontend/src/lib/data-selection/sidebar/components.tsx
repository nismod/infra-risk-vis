import { ArrowRight } from '@mui/icons-material';
import { Stack } from '@mui/material';
import { FC, createContext, useContext } from 'react';
import { useRecoilState } from 'recoil';

import { RecoilStateFamily } from '@/lib/recoil/types';

import { Accordion, AccordionDetails, AccordionSummary, AccordionTitle, ExpandablePanel } from './ExpandablePanel';
import { VisibilityToggle } from './VisibilityToggle';
import { SubPath, usePath } from './paths';

const VisibilityStateContext = createContext<RecoilStateFamily<boolean, string>>(null);

export function useVisibilityState(path: string) {
  const visibilityState = useContext(VisibilityStateContext);

  return useRecoilState(visibilityState(path));
}

const ExpandedStateContext = createContext<RecoilStateFamily<boolean, string>>(null);

export function useExpandedState(path: string) {
  const expandedState = useContext(ExpandedStateContext);
  return useRecoilState(expandedState(path));
}

export const Section: FC<{ path: string; title: string }> = ({ path, title, children }) => {
  return (
    <SubPath path={path}>
      <SectionImpl title={title}>{children}</SectionImpl>
    </SubPath>
  );
};

const SectionImpl = ({ title, children }) => {
  const path = usePath();
  const [visible, setVisible] = useVisibilityState(path);
  const [expanded, setExpanded] = useExpandedState(path);

  return (
    <Accordion
      title={title}
      expanded={expanded}
      onChange={(e, expanded) => setExpanded(expanded)}
      disableGutters
      sx={{
        bgcolor: '#fafafa',
        paddingLeft: '0px',
      }}
    >
      <AccordionSummary
        sx={(theme) => ({
          '& + .MuiCollapse-root': {
            borderLeft: '4px solid #fafafa',
          },
          '&:hover + .MuiCollapse-root': {
            borderLeftColor: theme.palette.primary.main,
          },
        })}
      >
        <AccordionTitle
          title={title}
          actions={
            <VisibilityToggle
              visibility={visible}
              onVisibility={(visible) => {
                setVisible(visible);
                setExpanded(visible);
              }}
            />
          }
        />
      </AccordionSummary>
      <AccordionDetails sx={{ padding: '0.5em', paddingRight: 0 }}>
        <Stack spacing={0.5}>{children}</Stack>
      </AccordionDetails>
    </Accordion>
  );
};

export const Layer: FC<{ path: string; title: string; disabled?: boolean }> = ({
  path,
  title,
  disabled = false,
  children,
}) => {
  return (
    <SubPath path={path}>
      <LayerImpl title={title} disabled={disabled}>
        {children}
      </LayerImpl>
    </SubPath>
  );
};

const LayerImpl: FC<{ title: string; disabled?: boolean }> = ({ title, disabled = false, children }) => {
  const path = usePath();
  const [visible, setVisible] = useVisibilityState(path);
  const [expanded, setExpanded] = useExpandedState(path);

  const allowExpand = visible && children != null;

  return (
    <Accordion
      disabled={disabled}
      title={title}
      expanded={allowExpand && expanded}
      onChange={(e, expanded) => setExpanded(expanded)}
      disableGutters
      sx={{
        border: '2px solid #eee',
      }}
      elevation={0}
    >
      <AccordionSummary
        sx={{ cursor: allowExpand ? 'pointer' : 'default' }}
        expandIcon={<ArrowRight color={allowExpand ? 'action' : 'disabled'} />}
      >
        <AccordionTitle
          title={title}
          actions={
            disabled ? null : (
              <VisibilityToggle
                visibility={visible}
                onVisibility={(visible) => {
                  setVisible(visible);
                  setExpanded(visible);
                }}
              />
            )
          }
        />
      </AccordionSummary>
      <AccordionDetails
        sx={{
          padding: 2,
          bgcolor: '#f5f5f5',
          border: '4px solid white',
        }}
      >
        {children}
      </AccordionDetails>
    </Accordion>
  );
  /*
  return (
    <ExpandablePanel
      disabled={disabled}
      title={title}
      expanded={expanded}
      onExpanded={setExpanded}
      allowExpand={visible && children != null}
      actions={
        <VisibilityToggle
          visibility={visible}
          onVisibility={(visible) => {
            setVisible(visible);
            setExpanded(visible);
          }}
        ></VisibilityToggle>
      }
    >
      {children}
    </ExpandablePanel>
  );
  */
};

export const SidebarPanel: FC<{ path: string; title: string }> = ({ path, title, children }) => {
  const [expanded, setExpanded] = useExpandedState(path);
  return (
    <ExpandablePanel title={title} expanded={expanded} onExpanded={setExpanded} actions={null}>
      {children}
    </ExpandablePanel>
  );
};

export const SidebarRoot: FC<{
  visibilityState: RecoilStateFamily<boolean, string>;
  expandedState: RecoilStateFamily<boolean, string>;
}> = ({ visibilityState, expandedState, children }) => {
  return (
    <VisibilityStateContext.Provider value={visibilityState}>
      <ExpandedStateContext.Provider value={expandedState}>{children}</ExpandedStateContext.Provider>
    </VisibilityStateContext.Provider>
  );
};
