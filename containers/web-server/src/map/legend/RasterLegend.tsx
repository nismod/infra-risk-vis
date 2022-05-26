import { RASTER_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { ViewLayer } from 'lib/data-map/view-layers';
import { FC, useCallback } from 'react';
import { GradientLegend } from './GradientLegend';
import { useRasterColorMapValues } from './use-color-map-values';

export const RasterLegend: FC<{ viewLayer: ViewLayer }> = ({ viewLayer }) => {
  const {
    params: { hazardType },
  } = viewLayer;
  const { label, dataUnit } = HAZARDS_METADATA[hazardType];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const { error, loading, colorMapValues } = useRasterColorMapValues(scheme, range);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()} ${dataUnit}`, [dataUnit]);

  return (
    <GradientLegend
      label={label}
      range={range}
      colorMapValues={!(error || loading) ? colorMapValues : null}
      getValueLabel={getValueLabel}
    />
  );
};
