import { FC, createContext, useContext } from 'react';
import { useRecoilState } from 'recoil';

import { RecoilStateFamily } from '@/lib/recoil/types';

import { ExpandablePanel } from './ExpandablePanel';
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
    <ExpandablePanel
      title={title}
      expanded={expanded}
      onExpanded={setExpanded}
      actions={
        <VisibilityToggle
          visibility={visible}
          onVisibility={(visible) => {
            setVisible(visible);
            setExpanded(visible);
          }}
        />
      }
    >
      {children}
    </ExpandablePanel>
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
