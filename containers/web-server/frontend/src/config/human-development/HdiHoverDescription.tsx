import { Typography } from '@mui/material';
import { FC } from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { DataDescription } from '@/lib/ui/data-display/DataDescription';

import { HDI_REGION_LEVEL_METADATA } from './metadata';

export const HdiHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const regionNameField = HDI_REGION_LEVEL_METADATA[viewLayer.params.regionLevel].nameField;

  return (
    <>
      <Typography variant="body2">{feature.properties[regionNameField]}</Typography>
      <DataDescription
        viewLayer={viewLayer}
        feature={feature}
        colorMap={viewLayer.styleParams?.colorMap}
      />
    </>
  );
};
