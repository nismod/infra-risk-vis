import { FC, useMemo } from 'react';

import { RASTER_COLOR_MAPS } from '../../../config/color-maps';
import { VIEW_LAYERS } from '../../../config/view-layers';
import { LAYERS } from '../../../config/layers';

import { RasterHover } from '../../DataMap';

import { useRasterColorMapValues } from '../../legend/use-color-map-values';

function useRasterColorMapLookup(colorMapValues) {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }) => [color, value])),
    [colorMapValues],
  );
}

export const RasterHoverDescription: FC<{ hoveredObject: RasterHover }> = ({ hoveredObject }) => {
  const { color } = hoveredObject;

  const { logicalLayer, deckLayer } = hoveredObject;
  const { label, dataUnit } = LAYERS[logicalLayer];
  const {
    dataParams: { hazardType },
  } = VIEW_LAYERS[deckLayer];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const title = `${label} (${dataUnit})`;

  const { colorMapValues } = useRasterColorMapValues(scheme, range);
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = `rgb(${color[0]},${color[1]},${color[2]})`;
  return (
    <div>
      <span style={{ color: colorString }}>■</span>&nbsp;
      <strong>{title}</strong>
      {rasterValueLookup?.[colorString] && <span>: {rasterValueLookup[colorString].toFixed(1)}</span>}
    </div>
  );
};
