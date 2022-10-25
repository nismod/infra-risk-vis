import { Box, Divider, Paper, Stack } from '@mui/material';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { ColorMap, FormatConfig } from '@/lib/data-map/view-layers';

import { viewLayersFlatState } from '@/state/layers/view-layers-flat';

import { VectorLegend } from './VectorLegend';

export const MapLegend: FC<{}> = () => {
  const viewLayers = useRecoilValue(viewLayersFlatState);

  const rasterLegends = [];

  const vectorLegendConfigs: Record<
    string,
    {
      colorMap: ColorMap;
      formatConfig: FormatConfig;
    }
  > = {};

  viewLayers.forEach((viewLayer) => {
    if (viewLayer.spatialType === 'raster') {
      // for raster layers, use the renderLegend() method of the view layer
      const layerLegend = viewLayer.renderLegend?.();
      if (layerLegend) {
        rasterLegends.push(layerLegend);
      }
    } else {
      const { colorMap } = viewLayer.styleParams ?? {};

      if (colorMap) {
        /**
         * Construct a key for grouping legends of the same type,
         * to avoid displaying the same legend many times for multiple layers.
         * Currently this is based on fieldGroup-field pair,
         * but this could need reworking for future cases.
         */
        const colorMapKey = `${colorMap.fieldSpec.fieldGroup}-${colorMap.fieldSpec.field}`;

        // save the colorMap and formatConfig for first layer of each group
        if (vectorLegendConfigs[colorMapKey] == null) {
          const { legendDataFormatsFn, dataFormatsFn } = viewLayer;

          const formatFn = legendDataFormatsFn ?? dataFormatsFn;
          const formatConfig = formatFn(colorMap.fieldSpec);

          vectorLegendConfigs[colorMapKey] = { colorMap, formatConfig };
        }
      }
    }
  });

  return rasterLegends.length || Object.keys(vectorLegendConfigs).length ? (
    <Paper>
      <Box p={1} maxWidth={270}>
        <Stack gap={0.3} divider={<Divider />}>
          {rasterLegends}
          {Object.entries(vectorLegendConfigs).map(([legendKey, { colorMap, formatConfig }]) => (
            <VectorLegend key={legendKey} colorMap={colorMap} legendFormatConfig={formatConfig} />
          ))}
        </Stack>
      </Box>
    </Paper>
  ) : null;
};
