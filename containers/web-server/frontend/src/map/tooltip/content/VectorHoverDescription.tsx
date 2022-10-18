import { Typography } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';

import { NETWORKS_METADATA } from '@/config/networks/metadata';
import { DataItem } from '@/details/features/detail-components';
import { singleViewLayerParamsState } from '@/state/layers/view-layers-params';

import { DataDescription } from '../DataDescription';
import { ColorBox } from './ColorBox';

export const VectorHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const layerParams = useRecoilValue(singleViewLayerParamsState(viewLayer.id));
  const { styleParams } = layerParams;
  const { colorMap } = styleParams ?? {};

  const isDataMapped = colorMap != null;

  const { label: title, color = '#ccc' } = NETWORKS_METADATA[viewLayer.params.assetId];

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
