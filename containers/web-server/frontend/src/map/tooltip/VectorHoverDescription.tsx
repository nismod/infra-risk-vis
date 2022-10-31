import { Typography } from '@mui/material';
import { FC } from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { ColorBox } from '@/lib/ui/data-display/ColorBox';
import { DataDescription } from '@/lib/ui/data-display/DataDescription';
import { DataItem } from '@/lib/ui/data-display/DataItem';

export const VectorHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
  metadataLookup: Record<string, any>;
}> = ({ hoveredObject, metadataLookup }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const { styleParams = {} } = viewLayer;
  const { colorMap } = styleParams;

  const isDataMapped = colorMap != null;

  const { label: title, color = '#ccc' } = metadataLookup[viewLayer.params.assetId];

  return (
    <>
      <Typography variant="body2">
        <ColorBox color={color} empty={isDataMapped} />
        {title}
      </Typography>

      <DataItem label="ID" value={feature.properties.asset_id} />
      {colorMap && <DataDescription viewLayer={viewLayer} feature={feature} colorMap={colorMap} />}
    </>
  );
};
