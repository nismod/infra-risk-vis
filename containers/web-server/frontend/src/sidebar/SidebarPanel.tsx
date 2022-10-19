import { useRecoilState } from 'recoil';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { ExpandablePanel } from '@/lib/ui/ExpandablePanel';
import { VisibilityToggle } from '@/lib/ui/VisibilityToggle';

import { sectionVisibilityState, sidebarSectionExpandedState } from '@/state/sections';

export const SidebarPanel = ({ id, title, children }) => {
  const [expanded, setExpanded] = useRecoilState(sidebarSectionExpandedState(id));
  const [visibility, setVisibility] = useRecoilState(sectionVisibilityState(id));

  return (
    <ExpandablePanel
      expanded={expanded}
      onExpanded={setExpanded}
      title={title}
      actions={<VisibilityToggle visibility={visibility} onVisibility={setVisibility} />}
    >
      <ErrorBoundary message="There was a problem displaying this section.">{children}</ErrorBoundary>
    </ExpandablePanel>
  );
};
