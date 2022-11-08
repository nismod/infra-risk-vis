import { ArrowRight } from '@mui/icons-material';
import { Stack } from '@mui/material';
import { FC, Suspense, createContext, forwardRef, useContext, useEffect } from 'react';
import { Flipped } from 'react-flip-toolkit';
import { useRecoilState } from 'recoil';

import { RecoilStateFamily } from '@/lib/recoil/types';

import { Accordion, AccordionDetails, AccordionSummary, AccordionTitle, ExpandablePanel } from './ExpandablePanel';
import { VisibilityToggle } from './VisibilityToggle';
import { PathContext, getSubPath, usePath } from './paths';

export const VisibilityStateContext = createContext<RecoilStateFamily<boolean, string>>(null);

export function useVisibilityState(path: string) {
  const visibilityState = useContext(VisibilityStateContext);

  return useRecoilState(visibilityState(path));
}

export const ExpandedStateContext = createContext<RecoilStateFamily<boolean, string>>(null);

export function useExpandedState(path: string) {
  const expandedState = useContext(ExpandedStateContext);
  return useRecoilState(expandedState(path));
}

export const PathChildrenStateContext = createContext<RecoilStateFamily<string[], string>>(null);

export function usePathChildrenState(path: string) {
  const pathChildrenState = useContext(PathChildrenStateContext);
  return useRecoilState(pathChildrenState(path));
}

function addValueToArray(val: string) {
  return (arr: string[]) => Array.from(new Set([...arr, val]));
}
function removeValueFromArray(val: string) {
  return (arr: string[]) => [...arr.filter((x) => x !== val)];
}

function useRegisterChild(parentPath, subPath) {
  const [, setPathChildren] = usePathChildrenState(parentPath);

  useEffect(() => {
    setPathChildren(addValueToArray(subPath));

    return () => {
      setPathChildren(removeValueFromArray(subPath));
    };
  }, [subPath, setPathChildren]);
}

export const SubPath: FC<{ path: string }> = ({ path, children }) => {
  const parentPath = usePath();
  const subPath = getSubPath(parentPath, path);

  useRegisterChild(parentPath, path);
  return <PathContext.Provider value={subPath}>{children}</PathContext.Provider>;
};

interface SectionProps {
  title: string;
}
export const Section: FC<{ path: string } & SectionProps> = ({ path, ...otherProps }) => {
  return (
    <SubPath path={path}>
      <SectionImpl {...otherProps} />
    </SubPath>
  );
};

const SectionImpl: FC<SectionProps> = ({ title, children }) => {
  const path = usePath();
  const [visible, setVisible] = useVisibilityState(path);
  const [expanded, setExpanded] = useExpandedState(path);

  return (
    <Flipped flipId={path}>
      <div>
        <Accordion
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
      </div>
    </Flipped>
  );
};

interface LayerProps {
  title: string;
  disabled?: boolean;
  unmountOnHide?: boolean;
}
export const Layer: FC<{ path: string } & LayerProps> = ({ path, ...otherProps }) => {
  return (
    <SubPath path={path}>
      <Suspense fallback="Loading layer data...">
        <LayerImpl {...otherProps} />
      </Suspense>
    </SubPath>
  );
};

const LayerImpl: FC<LayerProps> = ({ title, disabled = false, unmountOnHide = false, children }) => {
  const path = usePath();
  const [visible, setVisible] = useVisibilityState(path);
  const [expanded, setExpanded] = useExpandedState(path);

  const allowExpand = visible && children != null;

  return (
    <Accordion
      disabled={disabled}
      expanded={allowExpand && expanded}
      onChange={(e, expanded) => setExpanded(expanded)}
      disableGutters
      sx={{
        border: '2px solid #eee',
      }}
      elevation={0}
      TransitionProps={{
        unmountOnExit: !visible && unmountOnHide,
      }}
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

export const SidebarPanel: FC<{ path: string; title: string }> = ({ path, ...otherProps }) => {
  return (
    <SubPath path={path}>
      <SidebarPanelImpl {...otherProps} />
    </SubPath>
  );
};

const SidebarPanelImpl: FC<{ title: string }> = ({ title, children }) => {
  const path = usePath();
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
  pathChildrenState: RecoilStateFamily<string[], string>;
}> = ({ visibilityState, expandedState, pathChildrenState, children }) => {
  return (
    <VisibilityStateContext.Provider value={visibilityState}>
      <ExpandedStateContext.Provider value={expandedState}>
        <PathChildrenStateContext.Provider value={pathChildrenState}>{children}</PathChildrenStateContext.Provider>
      </ExpandedStateContext.Provider>
    </VisibilityStateContext.Provider>
  );
};
