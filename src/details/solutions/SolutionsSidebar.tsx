import { SidePanel } from 'details/SidePanel';
import { selectionState } from 'lib/data-map/interactions/interaction-state';
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
      <SolutionsSidebarContent feature={feature} solutionType={viewLayer.id} />
    </SidePanel>
  );
};
