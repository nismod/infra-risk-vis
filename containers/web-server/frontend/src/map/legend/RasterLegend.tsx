import { FC } from 'react';

import { GradientLegend } from './GradientLegend';
import { useRasterColorMapValues } from './use-color-map-values';

export interface RasterColorMap {
  scheme: string;
  range: [number, number];
}

export const RasterLegend: FC<{
  label: string;
  description?: string;
  colorMap: RasterColorMap;
  getValueLabel: (x: any) => string;
}> = ({ label, description, colorMap: { scheme, range }, getValueLabel }) => {
  const { error, loading, colorMapValues } = useRasterColorMapValues(scheme, range);

  return (
    <GradientLegend
      label={label}
      description={description}
      range={range}
      colorMapValues={!(error || loading) ? colorMapValues : null}
      getValueLabel={getValueLabel}
    />
  );
};
