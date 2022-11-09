import { FC, ReactNode } from 'react';

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
  getValueLabel: (x: any) => ReactNode | string;
}> = ({ label, description, colorMap: { scheme, range }, getValueLabel }) => {
  const colorMapValues = useRasterColorMapValues(scheme, range);

  return (
    <GradientLegend
      label={label}
      description={description}
      range={range}
      colorMapValues={colorMapValues}
      getValueLabel={getValueLabel}
    />
  );
};
