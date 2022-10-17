import { Box, Divider, Paper, Stack } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { ColorMap, FormatConfig } from '@/lib/data-map/view-layers';

import { viewLayersFlatState } from '@/state/layers/view-layers-flat';
import { viewLayersParamsState } from '@/state/layers/view-layers-params';

import { RasterLegend } from './RasterLegend';
import { VectorLegend } from './VectorLegend';

export const MapLegend: FC<{}> = () => {
  const viewLayers = useRecoilValue(viewLayersFlatState);
  const viewLayersParams = useRecoilValue(viewLayersParamsState);

  const hazardViewLayers = [];

  let dataColorMaps: Record<
    string,
    {
      colorMap: ColorMap;
      formatConfig: FormatConfig;
    }
  > = {};

  viewLayers.forEach((viewLayer) => {
    if (viewLayer.spatialType === 'raster') {
      hazardViewLayers.push(viewLayer);
    } else {
      /**
       * get style params from the viewLayerParams mechanism
       * (old mechanism for styleParams used by asset layers),
       * or the style params set directly in the new layer
       * (new mechanism used for styleParams by NBS, drought, population etc)
       */
      const { colorMap } = viewLayersParams[viewLayer.id].styleParams ?? viewLayer.styleParams ?? {};

      if (colorMap) {
        /**
         * Construct a key for grouping legends of the same type,
         * to avoid displaying the same legend many times for multiple layers.
         * Currently this is based on fieldGroup-field pair,
         * but this could need reworking for future cases.
         */
        const colorMapKey = `${colorMap.fieldSpec.fieldGroup}-${colorMap.fieldSpec.field}`;

        // save the colorMap and formatConfig for first layer of each group
        if (dataColorMaps[colorMapKey] == null) {
          const { legendDataFormatsFn, dataFormatsFn } = viewLayer;

          const formatFn = legendDataFormatsFn ?? dataFormatsFn;
          const formatConfig = formatFn(colorMap.fieldSpec);

          dataColorMaps[colorMapKey] = { colorMap, formatConfig };
        }
      }
    }
  });

  return hazardViewLayers.length || Object.keys(dataColorMaps).length ? (
    <Paper>
      <Box p={1} maxWidth={270}>
        <Stack gap={0.3} divider={<Divider />}>
          {hazardViewLayers.map((viewLayer) => (
            <RasterLegend key={viewLayer.id} viewLayer={viewLayer} />
          ))}
          {Object.entries(dataColorMaps).map(([legendKey, { colorMap, formatConfig }]) => (
            <VectorLegend key={legendKey} colorMap={colorMap} legendFormatConfig={formatConfig} />
          ))}
        </Stack>
      </Box>
    </Paper>
  ) : null;
};
