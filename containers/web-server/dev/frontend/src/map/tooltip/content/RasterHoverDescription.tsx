import { FC, useMemo } from 'react';

import { InteractionTarget, RasterTarget } from 'lib/data-map/interactions/use-interactions';

import { RASTER_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';

import { useRasterColorMapValues } from '../../legend/use-color-map-values';
import { ColorBox } from './ColorBox';
import { Box } from '@mui/material';
import { DataItem } from 'details/features/detail-components';

function useRasterColorMapLookup(colorMapValues) {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }) => [color, value])),
    [colorMapValues],
  );
}

function formatHazardValue(color, value, dataUnit) {
  return (
    <>
      <ColorBox color={color} />
      {value == null ? '' : value.toFixed(1) + dataUnit}
    </>
  );
}

export const RasterHoverDescription: FC<{ hoveredObject: InteractionTarget<RasterTarget> }> = ({ hoveredObject }) => {
  const { color } = hoveredObject.target;

  const {
    viewLayer: {
      id,
      params: { hazardType },
    },
  } = hoveredObject;
  const { label, dataUnit } = HAZARDS_METADATA[id];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const title = `${label}`;

  const { colorMapValues } = useRasterColorMapValues(scheme, range);
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = `rgb(${color[0]},${color[1]},${color[2]})`;
  const value = rasterValueLookup?.[colorString];
  return (
    <Box>
      <DataItem label={title} value={formatHazardValue(colorString, value, dataUnit)} />
    </Box>
  );
};
