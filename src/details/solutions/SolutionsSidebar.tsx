import { SidePanel } from 'details/SidePanel';
import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';
import { SolutionsSidebarContent } from './SolutionsSidebarContent';

export const SolutionsSidebar: FC<{}> = () => {
  const featureSelection = useRecoilValue(selectionState('solutions'));

  if (!featureSelection) return null;

  const {
    target: { feature },
    viewLayer,
  } = featureSelection;

  return (
    <SidePanel>
      <ErrorBoundary message="There was a problem displaying these details.">
        <SolutionsSidebarContent feature={feature} solutionType={viewLayer.id} />
      </ErrorBoundary>
    </SidePanel>
  );
};
