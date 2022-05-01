import { Box, Typography } from '@mui/material';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { ColorMap, FieldSpec, ViewLayer } from 'lib/data-map/view-layers';
import _ from 'lodash';
import { FC, useMemo } from 'react';
import { useRecoilValue } from 'recoil';
import { singleViewLayerParamsState } from 'state/layers/view-layers-params';
import { ColorBox } from './ColorBox';

export const DataDescription: FC<{
  viewLayer: ViewLayer;
  feature: any;
  colorMap: ColorMap;
}> = ({ viewLayer, feature, colorMap: { fieldSpec: colorField, colorSpec } }) => {
  const accessor = useMemo(() => viewLayer.dataAccessFn?.(colorField), [viewLayer, colorField]);

  const value = accessor?.(feature);

  const colorFn = useMemo(() => colorMap(colorSpec), [colorSpec]);

  const color = colorFn(value);

  const { getDataLabel, getValueFormatted } = viewLayer.dataFormatsFn(colorField);

  const dataLabel = getDataLabel(colorField);
  const formattedValue = getValueFormatted(value, colorField);

  return (
    <Box>
      <DataItem
        label={dataLabel}
        value={
          <>
            <ColorBox color={color} />
            {formattedValue ?? '-'}
          </>
        }
      />
    </Box>
  );
};

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
