import { FC, ReactNode } from 'react';
import { useRecoilState } from 'recoil';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { ExpandablePanel } from '@/lib/ui/ExpandablePanel';
import { VisibilityToggle } from '@/lib/ui/VisibilityToggle';

import { sectionVisibilityState, sidebarSectionExpandedState } from '@/state/sections';

export interface SidebarPanelProps {
  id: string;
  title: string;
  disabled?: boolean;
  children?: ReactNode;
}

export const SidebarPanel: FC<SidebarPanelProps> = ({ id, title, disabled = false, children = null }) => {
  const [expanded, setExpanded] = useRecoilState(sidebarSectionExpandedState(id));
  const [visibility, setVisibility] = useRecoilState(sectionVisibilityState(id));

  return (
    <ExpandablePanel
      expanded={expanded}
      onExpanded={setExpanded}
      disabled={disabled}
      title={title}
      actions={disabled ? null : <VisibilityToggle visibility={visibility} onVisibility={setVisibility} />}
    >
      <ErrorBoundary message="There was a problem displaying this section.">{children}</ErrorBoundary>
    </ExpandablePanel>
  );
};
