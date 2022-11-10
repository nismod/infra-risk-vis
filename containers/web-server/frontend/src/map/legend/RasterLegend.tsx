import { FC, ReactNode } from 'react';

import { GradientLegend } from './GradientLegend';
import { useRasterColorMapValues } from './use-color-map-values';

export interface RasterColorMap {
  scheme: string;
  range: [number, number];

  /**
   * If specified, determines whether the UI should indicate to the user that
   * a value at the min/max end of the range also represents values
   * below/above that end of the range.
   */
  rangeTruncated?: [boolean, boolean];
}

export interface ColorValue {
  color: string;
  value: any;
}
export interface RasterColorMapValues {
  colorMapValues: ColorValue[];
  rangeTruncated: [boolean, boolean];
}

export const RasterLegend: FC<{
  label: string;
  description?: string;
  colorMap: RasterColorMap;
  getValueLabel: (x: any) => ReactNode | string;
}> = ({ label, description, colorMap: { scheme, range, rangeTruncated = [false, false] }, getValueLabel }) => {
  const colorMapValues = useRasterColorMapValues(scheme, range);

  const colorMap = {
    colorMapValues,
    rangeTruncated,
  };

  return (
    <GradientLegend
      label={label}
      description={description}
      range={range}
      colorMap={colorMap}
      getValueLabel={getValueLabel}
    />
  );
};
