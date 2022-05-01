import React, { FC } from 'react';

import { FeatureSidebarContent } from './FeatureSidebarContent';
import { useRecoilValue } from 'recoil';
import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { SidePanel } from 'details/SidePanel';

export const FeatureSidebar: FC<{}> = () => {
  const featureSelection = useRecoilValue(selectionState('assets'));

  if (!featureSelection) return null;

  const {
    target: { feature },
    viewLayer,
  } = featureSelection;

  return (
    <SidePanel>
      <FeatureSidebarContent feature={feature} assetType={viewLayer.id} />
    </SidePanel>
  );
};
